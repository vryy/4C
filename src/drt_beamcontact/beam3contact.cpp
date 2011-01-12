/*!-----------------------------------------------------------------------------------------------------------
\file beam3contact.cpp
\brief One beam contact pair (two beam elements)

<pre>
Maintainer: Alexander Popp, Christian Cyron
            {popp,cyron}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*-----------------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "beam3contact.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#ifdef D_BEAM3
#include "../drt_beam3/beam3.H"
#endif
#ifdef D_BEAM3II
#include "../drt_beam3ii/beam3ii.H"
#endif

/*----------------------------------------------------------------------*
 |  constructor (public)                                      popp 04/10|
 *----------------------------------------------------------------------*/
CONTACT::Beam3contact::Beam3contact(const DRT::Discretization& pdiscret,
		                                const DRT::Discretization& cdiscret,
		                                const int& dofoffset,
		                                DRT::Element* element1,
                                    DRT::Element* element2,
                                    const Epetra_SerialDenseMatrix ele1pos,
                                    const Epetra_SerialDenseMatrix ele2pos):
pdiscret_(pdiscret),
cdiscret_(cdiscret),
dofoffset_(dofoffset),
element1_(element1),
element2_(element2),
ele1pos_(ele1pos),
ele2pos_(ele2pos),
contactflag_(false)
{
	// initialize augmented lagrange multiplier to zero
	lmuzawa_ = 0.0;
	
	// initialize class variables for contact point coordinates
	x1_.Size(NDIM);
	x2_.Size(NDIM);
	for(int i=0;i<NDIM;i++)
	{
		x1_[i]=0.0;
		x2_[i]=0.0;
	}
	
	return;
}

/*----------------------------------------------------------------------*
 |  copy-constructor (public)                                 popp 04/10|
 *----------------------------------------------------------------------*/
CONTACT::Beam3contact::Beam3contact(const Beam3contact& old):
pdiscret_(old.pdiscret_),
cdiscret_(old.cdiscret_),
dofoffset_(old.dofoffset_),
element1_(old.element1_),
element2_(old.element2_),
ele1pos_(old.ele1pos_),
ele2pos_(old.ele2pos_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate the element (public)                             popp 04/10|
 *----------------------------------------------------------------------*/
bool CONTACT::Beam3contact::Evaluate(LINALG::SparseMatrix& stiffmatrix,
                                     Epetra_Vector& fint, double& pp)
{
  //**********************************************************************
  // Evaluation of contact forces and stiffness
  //**********************************************************************
  // (1) Closest Point Projection
  //     -> handled via a local newton iteration
  //     -> results in the element coordinates xi and eta of contact point
	//     -> invalid porojections of pairs are sorted out
  // (2) Compute some auxiliary quantities
  //     -> normal vector, gap, shape functions, contact flag,
  //     -> linearizations of all geometric quantities
	// (3) Compute contact forces and stiffness
  //     -> stiffness terms are directly assembled to global matrix
  //     -> contact forces are only returned as global vector
  // (4) Perform some finite difference checks
  //     -> only if the flag BEAMCONTACTFDCHECKS is defined
	
  //**********************************************************************
  // (1) Closest Point Projection (CPP)
  //**********************************************************************
	// initialize the two element coordinates xi and eta of contact point
	vector<double> XiContact(2);
	
	// detect whether two elements lie on the same axis
	bool elementscolinear = false;
	
	// number of nodes of each element
	const int numnode1 = element1_->NumNode();
	const int numnode2 = element2_->NumNode();

	// vectors for shape functions and their derivatives
	Epetra_SerialDenseVector funct1(numnode1);					// = N1
	Epetra_SerialDenseVector funct2(numnode2);					// = N2
	Epetra_SerialDenseMatrix deriv1(1,numnode1);				// = N1,xi
	Epetra_SerialDenseMatrix deriv2(1,numnode2);				// = N2,eta
	Epetra_SerialDenseMatrix secondderiv1(1,numnode1);	// = N1,xixi
	Epetra_SerialDenseMatrix secondderiv2(1,numnode2);	// = N2,etaeta
	
	// coords and derivatives of the two contact points
	vector<double> x1(NDIM);														// = x1
	vector<double> x2(NDIM);														// = x2
	vector<double> dx1(NDIM);														// = x1,xi
	vector<double> dx2(NDIM);														// = x2,eta
	vector<double> ddx1(NDIM);													// = x1,xixi
	vector<double> ddx2(NDIM);													// = x2,etaeta

	// call function for closest point projection
	ClosestPointProjection(XiContact,elementscolinear);

	// check if the CPP found for this contact pair really valid
	// if not we are already done here...
	if (abs(XiContact[0])>1 || abs(XiContact[1])>1 || elementscolinear == true)
	{
		contactflag_ = false;
		return false;		
	}
	
	// debug output
	//cout << "\nContact between elements " << Element1().Id() << " and " << Element2().Id() << endl;
	//cout << "XiContact[0] = " << XiContact[0] << "	XiContact[1] = " << XiContact[1] << endl;
	
  //**********************************************************************
  // (2) Compute some auxiliary quantities
  //**********************************************************************
	// initialize
	vector<double> normal(3);
	double gap= 0.0;
	double norm = 0.0;
	
	// call function to fill variables for shape functions and their derivatives
	GetShapeFunctions(funct1,funct2,deriv1,deriv2,secondderiv1,secondderiv2,XiContact);
	
	// call function to fill variables with coords and derivs of the contact point
	ComputeCoordsAndDerivs(x1,x2,dx1,dx2,ddx1,ddx2,funct1,funct2,deriv1,deriv2,
	                       secondderiv1,secondderiv2,numnode1,numnode2);
	
	// call function to compute scaled normal and gap in possible contact point
	ComputeNormal(normal,gap,norm,x1,x2);
	
	// call function to evaluate contact status
	CheckContactStatus(pp);
	
  // store coordinates of contact point into class variables
  for(int i=0;i<NDIM;i++)
  {
    x1_[i]=x1[i];
    x2_[i]=x2[i];
  }
	
	//**********************************************************************
	// (3) Compute contact forces and stiffness
	//**********************************************************************
	
	// call function to evaluate and assemble contact forces
	EvaluateFcContact(pp,gap,normal,fint,funct1,funct2,numnode1,numnode2);
	
	// call function to evaluate and assemble contact stiffness
	EvaluateStiffcContact(pp,gap,normal,norm,stiffmatrix,x1,x2,dx1,dx2,
	                      ddx1,ddx2,funct1,funct2,deriv1,deriv2,
	                      secondderiv1,secondderiv2,numnode1,numnode2,XiContact);
	
  return true;
}

/*----------------------------------------------------------------------*
 |	Closest point projection                                  popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::ClosestPointProjection(vector<double>& XiContact, 
                                                   bool& elementscolinear)
{
	// local variables for element coordinates
	vector<double> eta(2);
		
	// number of nodes of each element
	const int numnode1 = element1_->NumNode();
	const int numnode2 = element2_->NumNode();
	
	// vectors for shape functions and their derivatives
	Epetra_SerialDenseVector funct1(numnode1);					// = N1
	Epetra_SerialDenseVector funct2(numnode2);					// = N2
	Epetra_SerialDenseMatrix deriv1(1,numnode1);				// = N1,xi
	Epetra_SerialDenseMatrix deriv2(1,numnode2);				// = N2,eta
	Epetra_SerialDenseMatrix secondderiv1(1,numnode1);	// = N1,xixi
	Epetra_SerialDenseMatrix secondderiv2(1,numnode2);	// = N2,etaeta
	
	// coords and derivatives of the two contacting points
	vector<double> x1(NDIM);														// = x1
	vector<double> x2(NDIM);														// = x2
	vector<double> dx1(NDIM);														// = x1,xi
	vector<double> dx2(NDIM);														// = x2,eta
	vector<double> ddx1(NDIM);													// = x1,xixi
	vector<double> ddx2(NDIM);													// = x2,etaeta
	
	// initialize function f and Jacobian df for Newton iteration
	vector<double> f(2);
	LINALG::Matrix<2,2> df;
	
	// initial scalar residual (L2-norm of f)
	double residual = 0.0;
	
	//**********************************************************************
	// local Newton iteration
	//**********************************************************************
	for (int i=0;i<BEAMCONTACTMAXITER;++i)			 
	{
		// reset shape function variables to zero
		for (int j=0;j<numnode1;j++) funct1[j]=0;
		for (int j=0;j<numnode2;j++) funct2[j]=0;
		for (int j=0;j<numnode1;j++) deriv1(0,j)=0;
		for (int j=0;j<numnode2;j++) deriv2(0,j)=0;
		for (int j=0;j<numnode1;j++) secondderiv1(0,j)=0;
		for (int j=0;j<numnode2;j++) secondderiv2(0,j)=0;
				
		// update shape functions and their derivatives
		GetShapeFunctions(funct1,funct2,deriv1,deriv2,secondderiv1,secondderiv2,eta);
		
		// update coordinates and derivatives of contact points
		ComputeCoordsAndDerivs(x1,x2,dx1,dx2,ddx1,ddx2,funct1,funct2,deriv1,deriv2,
		                       secondderiv1,secondderiv2,numnode1,numnode2);
		
		// compute norm of difference vector to scale the equations
		// (this yields better conditioning)
		double norm = sqrt((x1[0]-x2[0])*(x1[0]-x2[0])
		                 + (x1[1]-x2[1])*(x1[1]-x2[1])
		                 + (x1[2]-x2[2])*(x1[2]-x2[2]));
		
		// the closer the beams get, the smaller is norm
		// norm is not allowed to be too small, else numerical problems occur
		if (norm < 1e-12)	norm=1.0;
		
		// evaluate f at current eta[i]
		EvaluateNewtonF(f,x1,x2,dx1,dx2,norm);
		
		// compute the scalar redisuum
		residual = sqrt(f[0]*f[0] + f[1]*f[1]);
		
	  // check if Newton iteration has converged
		if (residual < BEAMCONTACTTOL) break;	
		
		// evaluate Jacobian of f at current eta[i]
		EvaluateNewtonGradF(df,x1,x2,dx1,dx2,ddx1,ddx2,norm);								
			
		// Inverting (2x2) matrix df' by hard coded formula, so that it is
		// possible to handle colinear vectors, because they lead to det=0
		double det_df = df(0,0)*df(1,1)-df(1,0)*df(0,1);
				
		//********************************************************************
		// ASSUMPTION:
		// If det_df=0 we assume, that the two elements have an identical
		// neutral axis. These contact objects will be rejected. The outcome
		// of this physically rare phenomenon is that handling of line contact
		// is not possible with this approach.
		//********************************************************************
		// NOTE:
		// For efficiency reasons it would be better to detect the colinearity
		// before creating the the contact object. This is not yet implemented.
		//********************************************************************
		// singular df
		if (det_df==0)
		{
		  // sort out
		  elementscolinear = true;
		}
		
		// regular df (inversion possible)
		else
		{
		  // do not sort out
		  elementscolinear = false;
		  
			// invert df
			LINALG::Matrix<2,2> df_inv;
			df_inv(0,0)=df(1,1)/det_df;
			df_inv(0,1)=-df(0,1)/det_df;
			df_inv(1,0)=-df(1,0)/det_df;
			df_inv(1,1)=df(0,0)/det_df;
		
			// update element coordinates of contact point
	    eta[0] += -df_inv(0,0)*f[0] - df_inv(0,1)*f[1];
	    eta[1] += -df_inv(1,0)*f[0] - df_inv(1,1)*f[1];
		}
	}
	//**********************************************************************
	
	// Newton iteration unconverged
  if (residual > BEAMCONTACTTOL && !elementscolinear)
  {
  	eta[0] = 1e+12;
  	eta[1] = 1e+12;
  }
  	  	
  // store and return final result
  XiContact[0]=eta[0];
  XiContact[1]=eta[1];
  
  return;	
}

/*----------------------------------------------------------------------*
 |	Evaluate function f in CPP                                popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::EvaluateNewtonF(vector<double>& f, const vector<double>& x1, 
		 const vector<double>& x2,	const vector<double>& dx1, const vector<double>& dx2, const double& norm)
{
	// reset f
	f[0]=0;
	f[1]=0;
	
	// evaluate f
	// see Wriggers, Computational Contact Mechanics, equation (12.5)
	for (int i=0;i<NDIM;i++)
	{
		f[0] += (x1[i]-x2[i])*dx1[i] / norm;
		f[1] += (x1[i]-x2[i])*dx2[i] / norm;
	}
	
	return;
}

/*----------------------------------------------------------------------*
 |	Evaluate Jacobian df in CPP                               popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::EvaluateNewtonGradF(LINALG::Matrix<2,2>& df, const vector<double>& x1,
		 const vector<double>& x2, const vector<double>& dx1, const vector<double>& dx2,
		 const vector<double>& ddx1, const vector<double>& ddx2, const double& norm)
{
	// reset df
	df(0,0)=0;
	df(0,1)=0;
	df(1,0)=0;
	df(1,1)=0;
	
	// evaluate df
	// see Wriggers, Computational Contact Mechanics, equation (12.7)
	for(int i=0;i<NDIM;i++)
	{
		df(0,0) += (dx1[i]*dx1[i] + (x1[i]-x2[i])*ddx1[i]) / norm;
		df(0,1) += dx1[i]*dx2[i] / norm;
		df(1,0) += -dx2[i]*dx1[i] / norm;
		df(1,1) += (-dx2[i]*dx2[i] + (x1[i]-x2[i])*ddx2[i]) / norm;
	}
	
	return;	
}

/*----------------------------------------------------------------------*
 |	Check if conact is active or inactive                     popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::CheckContactStatus(double& pp)
{
	// check contact condition valid for both pure penalty
  // (lmuzawa = 0) and augmented lagrange (lmuzawa != 0)
	if (lmuzawa_ - pp * gap_ > 0) contactflag_ = true;
	else                          contactflag_ = false;

	return;
}

/*----------------------------------------------------------------------*
 |	Get global dofs of a node                                 popp 12/10|
 *----------------------------------------------------------------------*/
vector<int> CONTACT::Beam3contact::GetGlobalDofs(DRT::Node* node)
{
  // get dofs in beam contact discretization
	vector<int> cdofs = ContactDiscret().Dof(node);

	// get dofs in problem discretization via offset
	vector<int> pdofs((int)(cdofs.size()));
	for (int k=0;k<(int)(cdofs.size());++k)
		pdofs[k] = cdofs[k]-DofOffset();

	return pdofs;
}

/*----------------------------------------------------------------------*
 |  Compute contact forces                                    popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::EvaluateFcContact(const double& pp,
     const double& gap, const vector<double>& normal, Epetra_Vector& fint,
     Epetra_SerialDenseVector funct1, Epetra_SerialDenseVector funct2,
		 const int numnode1, const int numnode2)	
{
	// get dimensions for vectors fc1 and fc2
	const int dim1 = NDIM*numnode1;
	const int dim2 = NDIM*numnode2;
	
	// temporary vectors for contact forces, DOF-GIDs and owning procs
	Epetra_SerialDenseVector fc1(dim1);
	Epetra_SerialDenseVector fc2(dim2);
	vector<int> lm1(dim1);
	vector<int> lm2(dim2);
	vector<int> lmowner1(dim1);
	vector<int> lmowner2(dim2);
	
	// flag indicating assembly
	bool DoNotAssemble = false;
	
	//**********************************************************************
	// evaluate contact forces for active pairs
	//**********************************************************************
	if (contactflag_)
	{
	  // node ids of both elements
		const int* node_ids1 = element1_->NodeIds();
		const int* node_ids2 = element2_->NodeIds();
		
		//********************************************************************
		// Compute Fc1 (force acting on first element)
		//********************************************************************
		for (int i=0;i<numnode1;++i)
		{
			// get node pointer and dof ids
			DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
			vector<int> NodeDofGIDs = GetGlobalDofs(node);
			
			// compute force vector Fc1 and prepare assembly
			for (int j=0;j<NDIM;++j)
			{
				fc1[NDIM*i+j] = (lmuzawa_ - pp*gap) * normal[j] * funct1[i];					
				lm1[NDIM*i+j] = NodeDofGIDs[j];
				lmowner1[NDIM*i+j] = node->Owner();
			}
		}
		
		//********************************************************************
		// Compute Fc2 (force acting on second element)
		//********************************************************************
		for(int i=0;i<numnode2;++i)
		{
			// get node pointer and dof ids
			DRT::Node* node = ContactDiscret().gNode(node_ids2[i]);
			vector<int> NodeDofGIDs =  GetGlobalDofs(node);
			
			// compute force vector Fc2 and prepare assembly
			for (int j=0;j<NDIM;++j)
			{
				fc2[NDIM*i+j] = -(lmuzawa_ - pp*gap) * normal[j] * funct2[i];	
				lm2[NDIM*i+j] = NodeDofGIDs[j];
				lmowner2[NDIM*i+j] = node->Owner();
			}
		}
	}
	
	//**********************************************************************
	// no forces for inactive pairs
	//**********************************************************************
	else
	{
	  // set flag to avoid assembly
		DoNotAssemble = true;
		
		// compute forces
		for (int i=0;i<NDIM*numnode1;++i) fc1[i]=0;
		for (int i=0;i<NDIM*numnode2;++i) fc2[i]=0;
	}
	
	//**********************************************************************
	// assemble contact forces
	//**********************************************************************
	if (!DoNotAssemble)
	{	
		// assemble fc1 and fc2 into global contact force vector
		LINALG::Assemble(fint,fc1,lm1,lmowner1);
  	LINALG::Assemble(fint,fc2,lm2,lmowner2);
  	
  	// debug output
  	//cout << "********************	FC *********************"<<endl;
		//for (int i=0;i<NDIM*numnode1;++i) cout << "Fc_1_" << i << ": " << fc1[i] << endl;	
		//for (int i=0;i<NDIM*numnode2;++i) cout << "Fc_2_" << i << ": " << fc2[i] << endl;	
		//cout << "******************* FINT ********************"<<endl;
		//for (int i=0;i<6*NDIM*(numnode1+numnode2);++i) cout << "fint_" << i << ": " << fint[i] << endl;
	}

	return;
}

/*----------------------------------------------------------------------*
 |	Evaluate contact stiffness                                popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::EvaluateStiffcContact(const double& pp,
     const double& gap, const vector<double>& normal, const double& norm, 
		 LINALG::SparseMatrix& stiffmatrix,const vector<double>& x1,const vector<double>& x2,
		 const vector<double>& dx1, const vector<double>& dx2,
		 const vector<double>& ddx1, const vector<double>& ddx2, 
		 const Epetra_SerialDenseVector& funct1, const Epetra_SerialDenseVector& funct2,	
		 const Epetra_SerialDenseMatrix& deriv1,	const Epetra_SerialDenseMatrix& deriv2,	
		 const Epetra_SerialDenseMatrix& secondderiv1,const Epetra_SerialDenseMatrix& secondderiv2, 
		 const int numnode1, const int numnode2, const vector<double>& XiContact)
{
	// temporary matrices for stiffness and vectors for DOF-GIDs and owning procs 
	Epetra_SerialDenseMatrix stiffc1(NDIM*numnode1,NDIM*(numnode1+numnode2));
	Epetra_SerialDenseMatrix stiffc2(NDIM*numnode2,NDIM*(numnode1+numnode2));
	vector<int> lmrow1(NDIM*numnode1);
	vector<int> lmrow2(NDIM*numnode2);
	vector<int> lmrowowner1(NDIM*numnode1);
	vector<int> lmrowowner2(NDIM*numnode2);
	vector<int> lmcol1(NDIM*(numnode1+numnode2));
	vector<int> lmcol2(NDIM*(numnode1+numnode2));
	
	// flag indicating assembly
	bool DoNotAssemble = false;
	
	//**********************************************************************
	// evaluate contact stiffness for active pairs
	//**********************************************************************
	if (contactflag_)
	{
		// auxiliary stiffmatrix for part III of linearization to avoid tensor notation
		Epetra_SerialDenseMatrix stiffc_III(NDIM*(numnode1+numnode2),NDIM*(numnode1+numnode2));
		
		// node ids of both elements
		const int* node_ids1 = element1_->NodeIds();
		const int* node_ids2 = element2_->NodeIds();
		
		// initialize storage for linearizations
		vector<double> delta_xi(NDIM*(numnode1+numnode2));
		vector<double> delta_eta(NDIM*(numnode1+numnode2));
		vector<double> distance(3);
		double normdist = 0.0;
		vector<double> delta_gap(NDIM*(numnode1+numnode2));
		Epetra_SerialDenseMatrix delta_x1_minus_x2(NDIM,NDIM*(numnode1+numnode2));
		Epetra_SerialDenseMatrix delta_n(NDIM, NDIM*(numnode1+numnode2));
		
		//********************************************************************
		// evaluate linearizations and distance
		//********************************************************************
		// linearization of contact point
		ComputeLinXiAndLinEta(delta_xi,delta_eta,x1,x2,dx1,dx2,ddx1,ddx2,funct1,funct2,
		                      deriv1,deriv2,normal,norm,numnode1,numnode2,XiContact);
		
		// evaluation of distance
		ComputeDistance(distance, normdist, normal, norm);
		
		// linearization of gap function
		ComputeLinGap(delta_gap,delta_xi,delta_eta,x1,x2,dx1,dx2,funct1,funct2,normdist,
		              numnode1,numnode2,normal,norm,gap,delta_x1_minus_x2);
		
		// linearization of normal vector
		ComputeLinNormal(delta_n,x1,x2,norm,numnode1,numnode2,delta_x1_minus_x2,normal,XiContact);
		
		//********************************************************************
		// prepare assembly
		//********************************************************************
		// fill lmrow1 and lmrowowner1
		for (int i=0;i<numnode1;++i)
		{
			// get pointer and dof ids
			DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
			vector<int> NodeDofGIDs =  GetGlobalDofs(node);
			
			for (int j=0;j<NDIM;++j)
			{
				lmrow1[NDIM*i+j]=NodeDofGIDs[j];
				lmrowowner1[NDIM*i+j]=node->Owner();
			}
		}
		
		// fill lmrow2 and lmrowowner2
		for (int i=0;i<numnode2;++i)
		{
			// get pointer and node ids
			DRT::Node* node = ContactDiscret().gNode(node_ids2[i]);
			vector<int> NodeDofGIDs =  GetGlobalDofs(node);
			
			for (int j=0;j<NDIM;++j)
			{
				lmrow2[NDIM*i+j]=NodeDofGIDs[j];
				lmrowowner2[NDIM*i+j]=node->Owner();
			}
		}
		
		// fill lmcol1 and lmcol2
		for (int i=0;i<numnode1;++i)
		{	
			// get pointer and node ids
			DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
			vector<int> NodeDofGIDs =  GetGlobalDofs(node);
			
			for (int j=0;j<NDIM;++j)
			{
				lmcol1[NDIM*i+j] = NodeDofGIDs[j];
				lmcol2[NDIM*i+j] = NodeDofGIDs[j];
			}
		}
		
		// fill lmcol1 and lmcol2
		for (int i=0;i<numnode2;++i)
		{
			// get pointer and node ids
			DRT::Node* node = ContactDiscret().gNode(node_ids2[i]);
			vector<int> NodeDofGIDs =  GetGlobalDofs(node);
			
			for (int j=0;j<NDIM;++j)
			{
			  lmcol1[NDIM*numnode1+NDIM*i+j] = NodeDofGIDs[j];
			  lmcol2[NDIM*numnode1+NDIM*i+j] = NodeDofGIDs[j];
			}
		}
				
		//********************************************************************
		// index vectors for access to shape function vectors and matrices
		//********************************************************************
		// In literature shape functions are stored in matrices but here only
		// a vector is used. To access an element of a shape function matrix,
		// use index vectors, that relate the vector index to the row index of
		// the shape function matrices. Use an if-contruction to select only
		// the entries of the shape function matrices, that are != 0.
		// --> if (j%NDIM == i): only quasi-diagonal entries are considered
		//********************************************************************
		// intialize index vectors
		vector<int> index1(NDIM*(numnode1+numnode2));
		vector<int> index2(NDIM*(numnode1+numnode2));
		
		// fill the index vectors
		for (int i=0;i<NDIM*numnode1;++i)
		{
			index1[i]=(int)floor(i/NDIM);
			index2[i]=(int)floor(i/NDIM);
		}
		for (int i=NDIM*numnode1;i<NDIM*(numnode1+numnode2);++i)
		{
			index1[i]=(int)floor((i-NDIM*numnode1)/NDIM);
			index2[i]=(int)floor((i-NDIM*numnode1)/NDIM);
		}
		
		//********************************************************************
		// evaluate contact stiffness
		// (1) stiffc1 of first element
		//********************************************************************
		
		//********************************************************************
		// part I
		//********************************************************************
		vector<double> normal_t_N1(NDIM*numnode1);
		for (int i=0;i<NDIM;i++)
		{
			for (int j=0;j<NDIM*numnode1;j++)
			{
				if (j%NDIM == i)
				{
					int row = index1[j];
					normal_t_N1[j] += normal[i] * funct1[row];
				}
			}
		}
		
		for (int i=0;i<NDIM*numnode1;i++)		
		{
			for (int j=0;j<NDIM*(numnode1+numnode2);j++)
			{
				stiffc1(i,j) = -pp * normal_t_N1[i] * delta_gap[j];
			}
		}
		
		//********************************************************************
		// part II
		//********************************************************************
		for  (int i=0;i<NDIM;i++)
		{
			for (int j=0;j<NDIM*numnode1;j++)
			{
				for (int k=0;k<NDIM*(numnode1+numnode2);k++)
				{
					if (j%NDIM == i)
					{
						int row = index1[j];
						stiffc1(j,k) += (lmuzawa_ - pp*gap) * funct1[row] * delta_n(i,k);
					}
				}
			}
		}

		//********************************************************************
		// part III
		//********************************************************************
		vector<double> normal_t_N1_xi(NDIM*numnode1);
		for (int i=0;i<NDIM;i++)
		{
			for (int j=0;j<NDIM*numnode1;j++)
			{
				if (j%NDIM == i)
				{
					int row = index1[j];
					normal_t_N1_xi[j] += normal[i] * deriv1(0,row);
				}
			}
		}
				
		for (int i=0;i<NDIM*numnode1;i++)
		{
			for (int j=0;j<NDIM*(numnode1+numnode2);j++)
			{
				stiffc1(i,j) += (lmuzawa_ - pp*gap) * normal_t_N1_xi[i] * delta_xi[j];
			}
		}
		
		//********************************************************************
	  // evaluate contact stiffness
		// (2) stiffc2 of second element
		//********************************************************************
		
		//********************************************************************
		// part I
		//********************************************************************
		vector<double> normal_t_N2(NDIM*numnode2);
		for (int i=0;i<NDIM;i++)
		{
			for (int j=0;j<NDIM*numnode2;j++)
			{
				if (j%NDIM == i)
				{
					int row = index2[j];
					normal_t_N2[j] += normal[i] * funct2[row];
				}
			}
		}
		
		for (int i=0;i<NDIM*numnode2;i++)		
		{
			for (int j=0;j<NDIM*(numnode1+numnode2);j++)
			{
				stiffc2(i,j) = pp * normal_t_N2[i] * delta_gap[j];
			}
		}
		
		//********************************************************************
		// part II
		//********************************************************************
		for (int i=0;i<NDIM;i++)
		{
			for (int j=0;j<NDIM*numnode2;j++)
			{
				for (int k=0;k<NDIM*(numnode1+numnode2);k++)
				{
					if (j%NDIM == i)
					{
						int row = index2[j];
						stiffc2(j,k) += -(lmuzawa_  - pp*gap) * funct2[row] * delta_n(i,k);
					}
				}
			}
		}

		//********************************************************************
		// part III
		//********************************************************************
		vector<double> normal_t_N2_eta(NDIM*numnode2);
		for (int i=0;i<NDIM;i++)
		{
			for (int j=0;j<NDIM*numnode2;j++)
			{
				if (j%NDIM == i)
				{
					int row = index2[j];
					normal_t_N2_eta[j] += normal[i] * deriv2(0,row);
				}
			}
		}

		for (int i=0;i<NDIM*numnode2;i++)
		{
			for (int j=0;j<NDIM*(numnode1+numnode2);j++)
			{
				stiffc2(i,j) += -(lmuzawa_ - pp*gap) * normal_t_N2_eta[i] * delta_eta[j];
			}
		}

		//********************************************************************
		// finite difference check of stiffc
		//********************************************************************
#ifdef BEAMCONTACTFDCHECKS		
		FDCheckStiffc(numnode1,numnode2,stiffc1,stiffc2,pp,normal,gap,funct1,funct2);
#endif
	}

	//**********************************************************************
  // no stiffness for inactive pairs
  //**********************************************************************
	else
	{
	   // set flag to avoid assembly
		DoNotAssemble = true;
		
    // compute stiffness
		for (int i=0;i<NDIM*numnode1;i++)
			for (int j=0;j<NDIM*numnode1;j++)
				stiffc1(i,j) = 0;
		for (int i=0;i<NDIM*numnode2;i++)
			for (int j=0;j<NDIM*numnode2;j++)
				stiffc2(i,j) = 0;
	}
	
	//**********************************************************************
  // assemble contact stiffness
  //**********************************************************************
	// change sign of stiffc1 and stiffc2 due to time integration.
  // according to analytical derivation there is no minus sign, but for
  // our time integration methods the negative stiffness must be assembled.
  for (int j=0;j<NDIM*(numnode1+numnode2);j++)
  {
    for (int i=0;i<NDIM*numnode1;i++)
      stiffc1(i,j) = stiffc1(i,j) * (-1);
    for (int i=0;i<NDIM*numnode2;i++)
      stiffc2(i,j) = stiffc2(i,j) * (-1);
  }
	  
	// now finally assemble stiffc1 and stiffc2
	if (!DoNotAssemble)
	{	
		stiffmatrix.Assemble(0,stiffc1,lmrow1,lmrowowner1,lmcol1);
		stiffmatrix.Assemble(0,stiffc2,lmrow2,lmrowowner2,lmcol2);
	}
	
	return;
}

/*----------------------------------------------------------------------*
 |	Compute normal vector in contact point                    popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::ComputeNormal(vector<double>& normal, double& gap,
    double& norm, const vector<double>& x1, const vector<double>& x2)
{
	// compute non-unit normal
	for (int i=0;i<NDIM;i++) normal[i] = x1[i]-x2[i];		

	// compute length of normal
	norm = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
	if (norm==0) dserror("ERROR: Normal of length zero!");
	
	// evlaute scalar gap function
	ComputeGap(gap,norm);
	
	// compute unit normal
	for (int i=0;i<NDIM;i++) normal[i] /= norm;
	
	return;
}

/*----------------------------------------------------------------------*
 |	Compute scalar gap funtion                                popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::ComputeGap(double& gap, const double& norm)
{
	// get moments of inertia of both elements
	// NOTE: here Iyy_ = Izz_ due to the circular cross section
	double MomentOfInertia_ele1 = 0;
	double MomentOfInertia_ele2 = 0;

  const DRT::ElementType & eot1 = element1_->ElementType();
  const DRT::ElementType & eot2 = element2_->ElementType();

  #ifdef D_BEAM3
      if ( eot1 == DRT::ELEMENTS::Beam3Type::Instance() )
        MomentOfInertia_ele1 = (static_cast<DRT::ELEMENTS::Beam3*>(element1_))->Iyy();
      if ( eot2 == DRT::ELEMENTS::Beam3Type::Instance() )
        MomentOfInertia_ele2 = (static_cast<DRT::ELEMENTS::Beam3*>(element2_))->Iyy();
  #endif  // #ifdef D_BEAM3
  #ifdef D_BEAM3II
      if ( eot1 == DRT::ELEMENTS::Beam3iiType::Instance() )
        MomentOfInertia_ele1 = (static_cast<DRT::ELEMENTS::Beam3ii*>(element1_))->Iyy();
      if ( eot2 == DRT::ELEMENTS::Beam3iiType::Instance() )
        MomentOfInertia_ele2 = (static_cast<DRT::ELEMENTS::Beam3ii*>(element2_))->Iyy();
  #endif  // #ifdef D_BEAM3II
	
	// compute radii of both elements
	double radius_ele1=0;
	double radius_ele2=0;
	ComputeEleRadius(radius_ele1,MomentOfInertia_ele1);
	ComputeEleRadius(radius_ele2,MomentOfInertia_ele2);
	

	// comute gap to be returned
	gap = norm - radius_ele1 - radius_ele2;
	
	// also set class variable
	gap_ = gap;

	return;
}

/*----------------------------------------------------------------------*
 |  compute radius from moment of inertia                     popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::ComputeEleRadius(double& radius, const double& moi)
{
  // fixed formula for circular cross sections
	radius = sqrt(sqrt(4 * moi / M_PI));
	
	return;
}

/*----------------------------------------------------------------------*
 | compute nodal coordinates and their derivatives            popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::ComputeCoordsAndDerivs(vector<double>& x1,
     vector<double>& x2, vector<double>& dx1, vector<double>& dx2,
     vector<double>& ddx1, vector<double>& ddx2,
     const Epetra_SerialDenseVector& funct1, const Epetra_SerialDenseVector& funct2,
		 const Epetra_SerialDenseMatrix& deriv1, const Epetra_SerialDenseMatrix& deriv2,
		 const Epetra_SerialDenseMatrix& secondderiv1, const Epetra_SerialDenseMatrix& secondderiv2,
		 const int& numnode1, const int& numnode2)
{
	// reset input variables
	for (int i=0;i<NDIM;i++)
	{
		x1[i]   = 0.0;
		x2[i]   = 0.0;
		dx1[i]  = 0.0;
		dx2[i]  = 0.0;
		ddx1[i] = 0.0;
		ddx2[i] = 0.0;
	}
		
	// auxialiary variables for the nodal coordinates
	Epetra_SerialDenseMatrix coord1(NDIM,numnode1);
	Epetra_SerialDenseMatrix coord2(NDIM,numnode2);
	
	// full coord1 and coord2
	for (int i=0;i<NDIM;i++)	
		for (int j=0;j<numnode1;j++)			
		  coord1(i,j)=ele1pos_(i,j);			
	for (int i=0;i<NDIM;i++)
		for (int j=0;j<numnode2;j++)	
		  coord2(i,j)=ele2pos_(i,j);
	
	// compute output variable
	for (int i=0;i<NDIM;i++)
	{
		for (int j=0;j<numnode1;j++)
			x1[i] += funct1[j] * coord1(i,j);							// x1 = N1 * x~1
		for (int j=0;j<numnode2;j++)
			x2[i] += funct2[j] * coord2(i,j);							// x2 = N2 * x~2
		for (int j=0;j<numnode1;j++)
			dx1[i] += deriv1(0,j) * coord1(i,j);					// dx1 = N1,xi * x~1
		for (int j=0;j<numnode2;j++)
			dx2[i] += deriv2(0,j) * coord2(i,j);					// dx2 = N2,eta * x~2
		for (int j=0;j<numnode1;j++)
			ddx1[i] += secondderiv1(0,j) * coord1(i,j);		// ddx1 = N1,xixi * x~1
		for (int j=0;j<numnode2;j++)
			ddx2[i] += secondderiv2(0,j) * coord2(i,j);		// ddx2 = N2,etaeta * x~2
	}
	
	return;
}

/*----------------------------------------------------------------------*
 |  evaluate shape functions and derivatives                  popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::GetShapeFunctions(
     Epetra_SerialDenseVector& funct1, Epetra_SerialDenseVector& funct2,
     Epetra_SerialDenseMatrix& deriv1, Epetra_SerialDenseMatrix& deriv2,
     Epetra_SerialDenseMatrix& secondderiv1, Epetra_SerialDenseMatrix& secondderiv2,
     const vector<double>& eta)
{
	// get both discretization types
	const DRT::Element::DiscretizationType distype1 = element1_->Shape();
	const DRT::Element::DiscretizationType distype2 = element2_->Shape();
			
	// get values and derivatives of shape functions
	DRT::UTILS::shape_function_1D(funct1, eta[0], distype1);
	DRT::UTILS::shape_function_1D(funct2, eta[1], distype2);
	DRT::UTILS::shape_function_1D_deriv1(deriv1, eta[0], distype1);
	DRT::UTILS::shape_function_1D_deriv1(deriv2, eta[1], distype2);
	DRT::UTILS::shape_function_1D_deriv2(secondderiv1, eta[0], distype1);
	DRT::UTILS::shape_function_1D_deriv2(secondderiv2, eta[1], distype2);
		
	return;
}

/*----------------------------------------------------------------------*
 |	Linearizations of contact point                           popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::ComputeLinXiAndLinEta(
    vector<double>& delta_xi, vector<double>& delta_eta, 
		const vector<double>& x1, const vector<double>& x2, const vector<double>& dx1, 
		const vector<double>& dx2, const vector<double>& ddx1, const vector<double>& ddx2,
		const Epetra_SerialDenseVector& funct1, const Epetra_SerialDenseVector& funct2, 
		const Epetra_SerialDenseMatrix& deriv1,	const Epetra_SerialDenseMatrix& deriv2, 
	  const vector<double>& normal, const double& norm,
	  const int numnode1, const int numnode2, const vector<double>& XiContact) 
{
  //**********************************************************************
  // we have to solve the following system of equations:
	//  _						  _		   _			_			 _							_			_			 _
	// | L(1,1)	L(1,2) |		| Lin_Xi  |		|	B(1,1)	B(1,2) |	 | Lin_d1 |
	// |							 | *  |					| = |								 | * |				|
	// |_L(2,1)	L(2,2)_|		|_Lin_Eta_|		|_B(2,1)	B(2,2)_|	 |_Lin_d2_|
	// 
	// this can be done easily because it is a linear 2x2-system.
  // we obtain the solution by inverting matrix L:
  // 
	// [Lin_Xi; Lin_Eta] = L^-1 * B * [Lin_d1; Lin_d2] = D * [Lin_d1; Lin_d2]
	// 
  //**********************************************************************
	
	// matrices to compute Lin_Xi and Lin_Eta
	Epetra_SerialDenseMatrix L(2,2);
	Epetra_SerialDenseMatrix L_inv(2,2);
	Epetra_SerialDenseMatrix B(2,NDIM*(numnode1+numnode2));
	Epetra_SerialDenseMatrix D(2,NDIM*(numnode1+numnode2));
	
	// compute L elementwise
	for (int i=0;i<NDIM;i++)
	{
		L(0,0) +=  dx1[i]*dx1[i] + (x1[i]-x2[i])*ddx1[i];
		L(0,1) += -dx2[i]*dx1[i];
		L(1,0) +=  dx1[i]*dx2[i];
		L(1,1) += -dx2[i]*dx2[i] + (x1[i]-x2[i])*ddx2[i];
	}
	
	// invert L by hand
	double det_L = L(0,0)*L(1,1) - L(0,1)*L(1,0);
	if (det_L == 0) dserror("ERROR: Determinant of L = 0");
	L_inv(0,0) =  L(1,1) / det_L;
	L_inv(0,1) = -L(0,1) / det_L;
	L_inv(1,0) = -L(1,0) / det_L;
	L_inv(1,1) =  L(0,0) / det_L;
	
	// index vectors for access to shape function vectors and matrices
	vector<int> index1(NDIM*(numnode1+numnode2));
	vector<int> index2(NDIM*(numnode1+numnode2));
	
	// fill the index vectors
	for (int i=0;i<NDIM*numnode1;i++)
	{
		index1[i]=(int) floor(i/3);
		index2[i]=(int) floor(i/3);
	}
	for (int i=NDIM*numnode1;i<NDIM*(numnode1+numnode2);i++)
	{
		index1[i]=(int) floor((i-NDIM*numnode1)/3);
		index2[i]=(int) floor((i-NDIM*numnode1)/3);
	}
	
	// compute B elementwise
	for (int i=0;i<NDIM*(numnode1+numnode2);i++)
	{
		int j=i%3;
		int row1 = index1[i];	
		int row2 = index2[i];
		
		// first block
		if (i<NDIM*numnode1)
		{
			B(0,i) = -(x1[j]-x2[j])*deriv1(0,row1) - dx1[j]*funct1[row1];
			B(1,i) = -dx2[j]*funct1[row1];
		}
		
		// second block
		else
		{
			B(0,i) = dx1[j]*funct2[row2];
			B(1,i) = -(x1[j]-x2[j])*deriv2(0,row2) + dx2[j]*funct2[row2];
		}
	}
			
	// compute D = L^-1 * B
	for (int i=0;i<2;i++)
		for (int j=0;j<2;j++)
			for (int k=0;k<NDIM*(numnode1+numnode2);k++)
				D(i,k) += L_inv(i,j) * B(j,k);		
	

	// finally the linearizations / directional derivatives
	for (int i=0;i<NDIM*(numnode1+numnode2);i++)
	{
		delta_xi[i] = D(0,i);
		delta_eta[i] = D(1,i);
	}

	// finite difference check
#ifdef BEAMCONTACTFDCHECKS
	FDCheckCPP(numnode1, numnode2, delta_xi, delta_eta,XiContact);
#endif
	
	return;
	
}

/*----------------------------------------------------------------------*
 |  Compute distance vector                                   popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::ComputeDistance(vector<double>& distance,
     double& normdist, const vector<double>& normal, const double& norm)
{
  // compute distance vector
	for (int i=0;i<(int)normal.size();i++)
		distance[i] = normal[i] * norm;
	
	// compute scalar distance
	normdist = norm;
	
	return;
}

/*----------------------------------------------------------------------*
 | Compute linearization of gap                               popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::ComputeLinGap(vector<double>& delta_gap,
      vector<double>& delta_xi, vector<double>& delta_eta,
      const vector<double>& x1, const vector<double>& x2, 
			const vector<double>& dx1, const vector<double>& dx2,
			const Epetra_SerialDenseVector& funct1, const Epetra_SerialDenseVector& funct2, 
			const double& normdist, const int& numnode1, const int& numnode2,
			const vector<double>& normal, const double& norm, const double& gap,
			Epetra_SerialDenseMatrix& delta_x1_minus_x2)
{
	// index vectors for access to shape function vectors and matrices
	vector<int> index1(NDIM*(numnode1+numnode2));
	vector<int> index2(NDIM*(numnode1+numnode2));
	
	// fill the index vectors
	for (int i=0;i<NDIM*numnode1;i++)
	{
		index1[i]=(int) floor(i/3);
		index2[i]=(int) floor(i/3);
	}
	
	for (int i=NDIM*numnode1;i<NDIM*(numnode1+numnode2);i++)
	{
		index1[i]=(int) floor((i-NDIM*numnode1)/3);
		index2[i]=(int) floor((i-NDIM*numnode1)/3);
	}
	
	// compute linearization of disctance vector
	for (int i=0;i<NDIM;i++)
	{
		for (int j=0;j<NDIM*(numnode1+numnode2);j++)
		{
		  // standard part for each j
			delta_x1_minus_x2(i,j) = dx1[i]*delta_xi[j] - dx2[i]*delta_eta[j];

			// only for first block
			if (j<NDIM*numnode1 && j%NDIM == i)
			{
				int row = index1[j];
				delta_x1_minus_x2(i,j) += funct1[row];
			}
			
			// only for second block
			else if(j%NDIM == i)
			{
				int row = index2[j];
				delta_x1_minus_x2(i,j) += -funct2[row];
			}
		}
	}
	
	// compute linearization of gap
	for (int i=0;i<NDIM;i++)								
		for (int j=0;j<NDIM*(numnode1+numnode2);j++)
			delta_gap[j] +=  (x1[i]-x2[i])/normdist * delta_x1_minus_x2(i,j);
	
	// finite difference check
#ifdef BEAMCONTACTFDCHECKS
  FDCheckLinGap(numnode1,numnode2,normal,norm,delta_gap,gap);
#endif
	
	return;	
}

/*----------------------------------------------------------------------*
 |  Compute linearization of normal                           popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::ComputeLinNormal(Epetra_SerialDenseMatrix& delta_n, 
		 const vector<double>& x1, const vector<double>& x2,
		 const double& norm,const int& numnode1, const int& numnode2,
		 const Epetra_SerialDenseMatrix& delta_x1_minus_x2, const vector<double>& normal, 
		 const vector<double>& XiContact)
{
	// local vectors for shape functions and their derivatives
	Epetra_SerialDenseVector funct1(numnode1);
	Epetra_SerialDenseVector funct2(numnode2);
	Epetra_SerialDenseMatrix deriv1(1,numnode1);
	Epetra_SerialDenseMatrix deriv2(1,numnode2);
	Epetra_SerialDenseMatrix secondderiv1(1,numnode1);
	Epetra_SerialDenseMatrix secondderiv2(1,numnode2);
	
	// get shape functions and their derivatives
	GetShapeFunctions(funct1,funct2,deriv1,deriv2,secondderiv1,secondderiv2,XiContact);
		
	// tensor product of normal (x) normal
	Epetra_SerialDenseMatrix n_tp_n(NDIM,NDIM);
	for (int i=0;i<NDIM;i++)
		for (int j=0;j<NDIM;j++)
			n_tp_n(i,j) = normal[i] * normal[j];

	// build a 3x3-identity matrix
	Epetra_SerialDenseMatrix identity(NDIM,NDIM);
	for (int i=0;i<NDIM;i++) identity(i,i) = 1;
			
	// compute linearization of normal
	for (int i=0;i<3;i++)
		for (int j=0;j<NDIM;j++)
			for (int k=0;k<NDIM*(numnode1+numnode2);k++)
				delta_n(i,k) += (identity(i,j) - n_tp_n(i,j)) * delta_x1_minus_x2(j,k) / norm;
	
	// finite difference check
#ifdef BEAMCONTACTFDCHECKS
  FDCheckLinNormal(numnode1,numnode2,delta_n,normal);
#endif	
	
	return;
}

/*----------------------------------------------------------------------*
 |  Reset Uzawa-based Lagrange mutliplier                   popp 04/2010|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::Resetlmuzawa()
{
	lmuzawa_ = 0.0;
	
	return;
}

/*----------------------------------------------------------------------*
 |  Update Uzawa-based Lagrange mutliplier                  popp 04/2010|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::Updatelmuzawa(const double& currentpp)
{
  // only update for active pairs, else reset
	if (contactflag_) lmuzawa_ -=  currentpp * GetGap();
	else              lmuzawa_ = 0.0;
	
	return;
}

/*----------------------------------------------------------------------*
 |  Update nodal coordinates (public)		 						          popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::UpdateElePos(Epetra_SerialDenseMatrix& newele1pos, 
		                                     Epetra_SerialDenseMatrix& newele2pos)
{
	ele1pos_ = newele1pos;
	ele2pos_ = newele2pos;
	
	return;
}

/*----------------------------------------------------------------------*
 |  FD check for closest point projection                     popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::FDCheckCPP(const int& numnode1, const int& numnode2, 
		 const vector<double>& delta_xi,	const vector<double>& delta_eta, 
		 const vector<double>& XiContact)
{
	// local boolean to detect whether two elements are colinear
	bool FD_elementscolinear = false;
	
	// local vectors for FD-approximations of delta_xi and delta_eta
	vector<double> FD_delta_xi(NDIM*(numnode1+numnode2));
	vector<double> FD_delta_eta(NDIM*(numnode1+numnode2));
	
	// step width and tolerances for FD-Check
	double eps = 1e-8;
	double rtol = 1e-2;					// relative tolerance
	double atol = 1e-3;					// absolute tolerance
	double toltreshold = 1e-6;	// analytical values < toltreshold are check with absolute tolerance
	
	// loop over ele1pos_
	for(int i=0;i<numnode1;i++)
	{
		for(int j=0;j<NDIM;j++)
		{
			// local vector with the two parameter coordinates xi1 and xi2 of the contact point
			vector<double> FD_XiContact;
			FD_XiContact.clear();
			FD_XiContact.resize(2);
				
			// step forward
			ele1pos_(j,i) += eps; 
			
			// do the closest point projection for the new set
			ClosestPointProjection(FD_XiContact,FD_elementscolinear);
			
			// FD-approximation
			FD_delta_xi[NDIM*i+j] = (FD_XiContact[0]-XiContact[0])/eps;
			FD_delta_eta[NDIM*i+j] = (FD_XiContact[1]-XiContact[1])/eps;
			
			// Undo step forward
			ele1pos_(j,i) -= eps;
		}
	}
	// loop over ele2pos_
	for(int i=0;i<numnode2;i++)
	{
		for(int j=0;j<NDIM;j++)
		{
			// local vector with the two parameter coordinates xi1 and xi2 of the contact point
			vector<double> FD_XiContact;
			FD_XiContact.clear();
			FD_XiContact.resize(2);
				
			// step forward
			ele2pos_(j,i) += eps; 
			
			// do the closest point projection for the new set
			ClosestPointProjection(FD_XiContact,FD_elementscolinear);
			
			// FD-approximation
			FD_delta_xi[NDIM*numnode1+NDIM*i+j] = (FD_XiContact[0]-XiContact[0])/eps;
			FD_delta_eta[NDIM*numnode1+NDIM*i+j] = (FD_XiContact[1]-XiContact[1])/eps;
			
			// Undo step forward
			ele2pos_(j,i) -= eps;
		}
	}
	
//	// Print the FD-Approximations to compare with analytical values
//	cout<<endl<<"Finite-Differences-Approximation of delta_xi and delta_eta:"<<endl<<endl;
//	for(int i=0;i<NDIM*(numnode1+numnode2);i++)
//	{					
//		cout << "FD_delta_xi_"<<i<<" = "<<FD_delta_xi[i]<<"	FD_delta_eta_"<<i<<" = "<<FD_delta_eta[i]<<endl;
//	}
	
	// Compare FD-Approximations and analytical values and give a warning, if deviation is to large
	for(unsigned int i=0;i<FD_delta_xi.size();i++)			// FD_delta_xi.size() = FD_delta_eta.size()
	{
		if(delta_xi[i] > toltreshold)		// then check relative error
		{
			if(fabs((FD_delta_xi[i]-delta_xi[i])/delta_xi[i]) > rtol)		// compare delta_xi elementwise
			{
				cout<<"FD_delta_xi_"<<i<<" = "<<FD_delta_xi[i]<<"	delta_xi_"<<i<<" = "<<delta_xi[i];
				cout<<"	********* WARNING from FD-Check: relative error > relative tolerance "<<rtol<<" *********"<<endl;
			}
		}
		else										// check absolute error
		{
			if(fabs(FD_delta_xi[i]-delta_xi[i]) > atol)		// compare delta_xi elementwise
			{
				cout<<"FD_delta_xi_"<<i<<" = "<<FD_delta_xi[i]<<"	delta_xi_"<<i<<" = "<<delta_xi[i];
				cout<<"	********* WARNING from FD-Check: absolute error "<<delta_xi[i]-FD_delta_xi[i];
				cout<<" > absolute tolerance "<<atol<<" *********"<<endl;
			}
		}
		if(delta_eta[i] > toltreshold)	// then check relative error
		{
			if(fabs((FD_delta_eta[i]-delta_eta[i])/delta_eta[i]) > rtol)		// compare delta_eta elementwise
			{
				cout<<"FD_delta_eta_"<<i<<" = "<<FD_delta_eta[i]<<"	delta_eta_"<<i<<" = "<<delta_eta[i];
				cout<<"	********* WARNING from FD-Check: relative error > relative tolerance "<<rtol<<" *********"<<endl;
			}
		}
		else										// check absolute error
		{
			if(fabs(FD_delta_eta[i]-delta_eta[i]) > atol)		// compare delta_eta elementwise
			{
				cout<<"FD_delta_eta_"<<i<<" = "<<FD_delta_eta[i]<<"	delta_eta_"<<i<<" = "<<delta_eta[i];
				cout<<"	********* WARNING from FD-Check: absolute error "<<delta_eta[i]-FD_delta_eta[i];
				cout<<" > absolute tolerance "<<atol<<" *********"<<endl;
			}
		}
	}
		
	return;
}

/*----------------------------------------------------------------------*
 |  FD check for gap function                                 popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::FDCheckLinGap(const int& numnode1, const int& numnode2, 
		 const vector<double>& normal, const double& norm,
		 const vector<double>& delta_gap, const double& gap)
{
	// local auxiliary variables for FD-approximations of delta_gap
	vector<double> FD_delta_gap(NDIM*(numnode1+numnode2));
	double FD_gap = 0.0;
	vector<double> FD_normal(3);
	double FD_norm = 0.0;
	vector<double> FD_XiContact(2);
	bool FD_elementscolinear = false;
	
	// step width and tolerances for FD-check
	double eps = 1e-10;
	double rtol = 1e-3;					// relative tolerance
	double atol = 1e-5;					// absolute tolerance
	double toltreshold = 1e-6;	// analytical values < toltreshold are check with absolute tolerance
		
	// local vectors for evaluated shape functions and their derivatives
	Epetra_SerialDenseVector FD_funct1(numnode1);
	Epetra_SerialDenseVector FD_funct2(numnode2);
	Epetra_SerialDenseMatrix FD_deriv1(1,numnode1);
	Epetra_SerialDenseMatrix FD_deriv2(1,numnode2);
	Epetra_SerialDenseMatrix FD_secondderiv1(1,numnode1);
	Epetra_SerialDenseMatrix FD_secondderiv2(1,numnode2);
		
	// local coords and derivatives of the two contacting points
	vector<double> FD_x1(NDIM);		// = x1
	vector<double> FD_x2(NDIM);		// = x2
	vector<double> FD_dx1(NDIM);		// = x1,xi
	vector<double> FD_dx2(NDIM);		// = x2,eta
	vector<double> FD_ddx1(NDIM);	// = x1,xixi
	vector<double> FD_ddx2(NDIM);	// = x2,etaeta
	
	// loop over ele1pos_
	for(int i=0;i<numnode1;i++)
	{
		for(int j=0;j<NDIM;j++)
		{
			// step forward
			ele1pos_(j,i) += eps; 
	
			// Get CPP-Solution for new set
			ClosestPointProjection(FD_XiContact, FD_elementscolinear);			
	
			// Call function to fill variables for shape functions and their derivatives
			GetShapeFunctions(FD_funct1,FD_funct2,FD_deriv1,FD_deriv2,FD_secondderiv1,FD_secondderiv2,FD_XiContact);
	
			// Call function to fill variabeles with coordinates and derivatives of the two contacting points
			ComputeCoordsAndDerivs(FD_x1,FD_x2,FD_dx1,FD_dx2,FD_ddx1,FD_ddx2,FD_funct1,FD_funct2,FD_deriv1,FD_deriv2,FD_secondderiv1,FD_secondderiv2,numnode1,numnode2);
	
			// Call function to compute scaled normal and gap in possible contact point
			ComputeNormal(FD_normal,FD_gap,FD_norm,FD_x1,FD_x2);
	
//			cout<<"FD_gap_"<<NDIM*i+j<<" = "<<FD_gap<<endl;
	
			// FD-approximation
			FD_delta_gap[NDIM*i+j] = (FD_gap - gap)/eps;
//			cout<<"FD_delta_gap_"<<NDIM*i+j<<" = "<<FD_delta_gap[NDIM*i+j]<<endl;
			
			// Undo step forward
			ele1pos_(j,i) -= eps;
		}
	}
	// loop over ele2pos_
	for(int i=0;i<numnode2;i++)
	{
		for(int j=0;j<NDIM;j++)
		{
			// step forward
			ele2pos_(j,i) += eps; 
			
			// Get CPP-Solution for new set
			ClosestPointProjection(FD_XiContact, FD_elementscolinear);			
		
			// Call function to fill variables for shape functions and their derivatives
			GetShapeFunctions(FD_funct1,FD_funct2,FD_deriv1,FD_deriv2,FD_secondderiv1,FD_secondderiv2,FD_XiContact);
			
			// Call function to fill variabeles with coordinates and derivatives of the two contacting points
			ComputeCoordsAndDerivs(FD_x1,FD_x2,FD_dx1,FD_dx2,FD_ddx1,FD_ddx2,FD_funct1,FD_funct2,FD_deriv1,FD_deriv2,FD_secondderiv1,FD_secondderiv2,numnode1,numnode2);
			
			// Call function to compute scaled normal and gap in possible contact point
			ComputeNormal(FD_normal,FD_gap,FD_norm,FD_x1,FD_x2);
			
//			cout<<"FD_gap_"<<NDIM*numnode1+NDIM*i+j<<" = "<<FD_gap<<endl;
			
			// FD-approximation
			FD_delta_gap[NDIM*numnode1+NDIM*i+j] = (FD_gap - gap)/eps;
//			cout<<"FD_delta_gap_"<<NDIM*numnode1+NDIM*i+j<<" = "<<FD_delta_gap[NDIM*numnode1+NDIM*i+j]<<endl;							
			
			// Undo step forward
			ele2pos_(j,i) -= eps;
		}
	}
	
//	// Print the FD-Approximations to compare with analytical values
//	cout<<endl<<"Finite-Differences-Approximation of delta_gap:"<<endl<<endl;
//	for(int i=0;i<NDIM*(numnode1+numnode2);i++)
//	{					
//		cout << "FD_delta_gap_"<<i<<" = "<<FD_delta_gap[i]<<endl;
//	}
		
	// Compare FD-Approximations and analytical values and give a warning, if deviation is to large
	cout<<endl;
	for(unsigned int i=0;i<FD_delta_gap.size();i++)
	{
		if(delta_gap[i] > toltreshold)	// then check relative error
		{
			if(fabs((FD_delta_gap[i]-delta_gap[i])/delta_gap[i]) > rtol)		// compare delta_gap elementwise
			{
				cout<<"FD_delta_gap_"<<i<<" = "<<FD_delta_gap[i]<<"	delta_gap_"<<i<<" = "<<delta_gap[i];
				cout<<"	********* WARNING from FD-Check: relative error > relative tolerance "<<rtol<<" *********"<<endl;
			}
		}	
		else														// check absolute error	
		{
			if(fabs(FD_delta_gap[i]-delta_gap[i]) > atol)
			{
				cout<<"FD_delta_gap_"<<i<<" = "<<FD_delta_gap[i]<<"	delta_gap_"<<i<<" = "<<delta_gap[i];
				cout<<"	********* WARNING from FD-Check: absolute error "<<FD_delta_gap[i]-delta_gap[i];
				cout<<" > absolute tolerance "<<atol<<" *********"<<endl;
			}
		}
	}
	
	return;
}

/*----------------------------------------------------------------------*
 |  FD check for normal vector                                popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::FDCheckLinNormal(const int& numnode1, const int& numnode2, 
		 const Epetra_SerialDenseMatrix& delta_n, const vector<double>& normal)
{
	// local auxiliary variables for FD-approximations of delta_n
	Epetra_SerialDenseMatrix FD_delta_n(NDIM, NDIM*(numnode1+numnode2));
		
	// step width and tolerances for FD-Check
	double eps = 1e-6;
	double rtol = 1e-2;					// relative tolerance
	double atol = 1e-3;					// absolute tolerance
	double toltreshold = 1e-6;	// analytical values < toltreshold are check with absolute tolerance
			
	// local vectors for evaluated shape functions and their derivatives
	Epetra_SerialDenseVector FD_funct1(numnode1);
	Epetra_SerialDenseVector FD_funct2(numnode2);
	Epetra_SerialDenseMatrix FD_deriv1(1,numnode1);
	Epetra_SerialDenseMatrix FD_deriv2(1,numnode2);
	Epetra_SerialDenseMatrix FD_secondderiv1(1,numnode1);
	Epetra_SerialDenseMatrix FD_secondderiv2(1,numnode2);
		
	// local coords and derivatives of the two contacting points
	vector<double> FD_x1(NDIM);		// = x1
	vector<double> FD_x2(NDIM);		// = x2
	vector<double> FD_dx1(NDIM);		// = x1,xi
	vector<double> FD_dx2(NDIM);		// = x2,eta
	vector<double> FD_ddx1(NDIM);	// = x1,xixi
	vector<double> FD_ddx2(NDIM);	// = x2,etaeta
	
//	cout<<endl<<"loop over ele1pos_"<<endl<<endl;
	
	// loop over ele1pos_
	for(int j=0;j<numnode1;j++)
	{
		for(int i=0;i<NDIM;i++)
		{
			vector<double> FD_normal(3);
			double FD_norm = 0.0;
			double FD_gap = 0.0;
			vector<double> FD_XiContact(2);
			bool FD_elementscolinear = false;
			
			// step forward
			ele1pos_(i,j) += eps;
//			cout<<"ele1pos_"<<ele1pos_<<endl;
	
			// Get CPP-Solution for new set
			ClosestPointProjection(FD_XiContact, FD_elementscolinear);			
//			cout<<"FD_XiContact:	"<<FD_XiContact[0]<<"	"<<FD_XiContact[1]<<endl;
			
			// Call function to fill variables for shape functions and their derivatives
			GetShapeFunctions(FD_funct1,FD_funct2,FD_deriv1,FD_deriv2,FD_secondderiv1,FD_secondderiv2,FD_XiContact);
	
			// Call function to fill variables with coordinates and derivatives of the two contacting points
			ComputeCoordsAndDerivs(FD_x1,FD_x2,FD_dx1,FD_dx2,FD_ddx1,FD_ddx2,FD_funct1,FD_funct2,FD_deriv1,FD_deriv2,FD_secondderiv1,FD_secondderiv2,numnode1,numnode2);
			
			// Test Print of FD_x1 and FD_x2
//			cout<<"FD_x1:	"<<FD_x1[0]<<"	"<<FD_x1[1]<<"	"<<FD_x1[2]<<endl;
//			cout<<"FD_x2:	"<<FD_x2[0]<<"	"<<FD_x2[1]<<"	"<<FD_x2[2]<<endl;
	
			// Call function to compute scaled normal and gap in possible contact point
			ComputeNormal(FD_normal,FD_gap,FD_norm,FD_x1,FD_x2);
						
//			// Test print of FD_normal
//			cout<<"FD_normal		"<<FD_normal[0]<<"	"<<FD_normal[1]<<"	"<<FD_normal[2]<<endl<<endl;
			
			for(int k=0;k<NDIM;k++)
			{
				// FD-approximation
				FD_delta_n(k,NDIM*j+i) += (FD_normal[k] - normal[k])/eps;
			}
			
//			cout<<"FD_delta_n_row"<<NDIM*j+i<<":	"<<FD_delta_n(0,NDIM*j+i)<<"	"<<FD_delta_n(1,NDIM*j+i)<<"	"<<FD_delta_n(2,NDIM*j+i)<<endl;
			
			// Undo step forward
			ele1pos_(i,j) -= eps;
		}
	}
	
//	cout<<endl<<"loop over ele2pos_"<<endl<<endl;
	
	// loop over ele2pos_
	for(int j=0;j<numnode2;j++)
	{
		for(int i=0;i<NDIM;i++)
		{
			vector<double> FD_normal(3);
			double FD_norm = 0.0;
			double FD_gap = 0.0;
			vector<double> FD_XiContact(2);
			bool FD_elementscolinear = false;
		
			// step forward
			ele2pos_(i,j) += eps; 
//			cout<<"ele2pos_"<<ele2pos_<<endl;
			
			// Get CPP-Solution for new set
			ClosestPointProjection(FD_XiContact, FD_elementscolinear);			
//			cout<<"FD_XiContact:	"<<FD_XiContact[0]<<"	"<<FD_XiContact[1]<<endl;
			
			// Call function to fill variables for shape functions and their derivatives
			GetShapeFunctions(FD_funct1,FD_funct2,FD_deriv1,FD_deriv2,FD_secondderiv1,FD_secondderiv2,FD_XiContact);
			
			// Call function to fill variables with coordinates and derivatives of the two contacting points
			ComputeCoordsAndDerivs(FD_x1,FD_x2,FD_dx1,FD_dx2,FD_ddx1,FD_ddx2,FD_funct1,FD_funct2,FD_deriv1,FD_deriv2,FD_secondderiv1,FD_secondderiv2,numnode1,numnode2);
			
//			for(int k=0;k<NDIM;k++)
//			{
//				FD_x1[k]=0;
//				FD_x2[k]=0;
//			}
//			for(int k=0;k<3;k++)		// Compute Coords of Contact Point after FD-manipulation manually
//			{
//				for(int m=0;m<numnode1;m++)
//					FD_x1[k] += FD_funct1[m] * ele1pos_(k,m);							// x1 = N1 * x~1			x~1 are the discret nodal coordinates
//				for(int m=0;m<numnode2;m++)
//					FD_x2[k] += FD_funct2[m] * ele2pos_(k,m);							// x2 = N2 * x~2
//			}
			
			// Test Print of FD_x1 and FD_x2
//			cout<<"FD_x1:	"<<FD_x1[0]<<"	"<<FD_x1[1]<<"	"<<FD_x1[2]<<endl;
//			cout<<"FD_x2:	"<<FD_x2[0]<<"	"<<FD_x2[1]<<"	"<<FD_x2[2]<<endl;
			
			// Call function to compute scaled normal and gap in possible contact point
			ComputeNormal(FD_normal,FD_gap,FD_norm,FD_x1,FD_x2);
			
//			// Test print of FD_normal
//			cout<<"FD_normal		"<<FD_normal[0]<<"	"<<FD_normal[1]<<"	"<<FD_normal[2]<<endl<<endl;
			
			for(int k=0;k<NDIM;k++)
			{
				// FD-approximation
				FD_delta_n(k,NDIM*numnode1+NDIM*j+i) += (FD_normal[k] - normal[k])/eps;
			}
			
//			cout<<"FD_delta_n_row"<<NDIM*numnode1+NDIM*j+i<<":	"<<FD_delta_n(0,NDIM*numnode1+NDIM*j+i)<<"	"<<FD_delta_n(1,NDIM*numnode1+NDIM*j+i)<<"	"<<FD_delta_n(2,NDIM*numnode1+NDIM*j+i)<<endl;
			
			// Undo step forward
			ele2pos_(i,j) -= eps;
		}
	}
	
//	cout<<"After FD-Approximation"<<endl;
	
//	// Print the FD-Approximations to compare with analytical values
//	cout<<endl<<"Finite-Differences-Approximation of delta_n:"<<endl<<endl;
//	cout<<"FD_delta_n:"<<endl<<FD_delta_n<<endl;
		
	// Compare FD-Approximations and analytical values and give a warning, if deviation is to large
	cout<<endl;
	for(int i=0;i<NDIM;i++)
	{
		for(int j=0;j<NDIM*(numnode1+numnode2);j++)
		{
			if(delta_n(i,j) > toltreshold)			// then check relative error
			{
				if(fabs((FD_delta_n(i,j)-delta_n(i,j))/delta_n(i,j)) > rtol)		// compare delta_n elementwise
				{
					cout<<"FD_delta_n_"<<i<<"_"<<j<<" = "<<FD_delta_n(i,j)<<"	delta_n_"<<i<<"_"<<j<<" = ";
					cout<<delta_n(i,j)<<"	********* WARNING from FD-Check: relative error = ";
					cout<<fabs((FD_delta_n(i,j)-delta_n(i,j))/delta_n(i,j))<<" > relative tolerance "<<rtol;
					cout<<" *********"<<endl;
				}
			}
			else																// check absolute error
			{
				if(fabs(FD_delta_n(i,j)-delta_n(i,j)) > atol)		// compare delta_n elementwise
				{
					cout<<"FD_delta_n_"<<i<<"_"<<j<<" = "<<FD_delta_n(i,j)<<"	delta_n_"<<i<<"_"<<j<<" = ";
					cout<<delta_n(i,j)<<"	********* WARNING from FD-Check: absolute error = ";
					cout<<fabs(FD_delta_n(i,j)-delta_n(i,j))<<" > absolute tolerance "<<atol;
					cout<<" *********"<<endl;
				}
			}
		}
	}
	
	return;
}

/*----------------------------------------------------------------------*
 |  FD check for contact stiffness                            popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3contact::FDCheckStiffc(const int& numnode1, const int& numnode2, 
     const Epetra_SerialDenseMatrix& stiffc1, const Epetra_SerialDenseMatrix& stiffc2,
     const double& pp, const vector<double>& normal, const double& gap,
     const Epetra_SerialDenseVector& funct1, const Epetra_SerialDenseVector& funct2)
{
	// local auxiliary variables for FD-approximations of stiffc1 and stiffc2
	vector<double> FD_delta_gap(NDIM*(numnode1+numnode2));
	double FD_gap = 0.0;
	vector<double> FD_normal(3);
	double FD_norm = 0.0;
	vector<double> FD_XiContact(2);
	bool FD_elementscolinear = false;
	
	// local matrices for FD-approximations of stiffc1 and stiffc2
	Epetra_SerialDenseMatrix FD_stiffc1(NDIM*numnode1,NDIM*(numnode1+numnode2));
	Epetra_SerialDenseMatrix FD_stiffc2(NDIM*numnode2,NDIM*(numnode1+numnode2));
	
	// step width and tolerances for FD-Check
	double eps = 1e-8;
	double rtol = 1e-3;					// relative tolerance
	double atol = 1e+2;					// absolute tolerance
	double toltreshold = 1e-6;	// analytical values < toltreshold are check with absolute tolerance
	
	// local vectors for evaluated shape functions and their derivatives
	Epetra_SerialDenseVector FD_funct1(numnode1);
	Epetra_SerialDenseVector FD_funct2(numnode2);
	Epetra_SerialDenseMatrix FD_deriv1(1,numnode1);
	Epetra_SerialDenseMatrix FD_deriv2(1,numnode2);
	Epetra_SerialDenseMatrix FD_secondderiv1(1,numnode1);
	Epetra_SerialDenseMatrix FD_secondderiv2(1,numnode2);
		
	// local coords and derivatives of the two contacting points
	vector<double> FD_x1(NDIM);		// = x1
	vector<double> FD_x2(NDIM);		// = x2
	vector<double> FD_dx1(NDIM);	// = x1,xi
	vector<double> FD_dx2(NDIM);	// = x2,eta
	vector<double> FD_ddx1(NDIM);	// = x1,xixi
	vector<double> FD_ddx2(NDIM);	// = x2,etaeta
	
	// Compute fc1 and fc2 analytical here explicit and not via the function "EvaluateFcContact"
	// because there may not be a second assembly of fc1 and fc2 
	
	// temporary vectors for contact forces, DOF-GIDs and owning procs 
	Epetra_SerialDenseVector fc1(NDIM*numnode1); // additional internal forces for element1 based on contact
	Epetra_SerialDenseVector fc2(NDIM*numnode2); // additional internal forces for element2 based on contact		
  
	for(int i=0;i<fc1.Length();i++)
		fc1[i]=0;
	for(int i=0;i<fc2.Length();i++)
		fc2[i]=0;
		
	if ((-gap)>0)		// Fc > 0 only, if penetration occurs
	{
		
		// Compute Fc1
		for(int i=0;i<numnode1;i++)
		{
			for(int j=0;j<NDIM;j++)
			{
				fc1[NDIM*i+j]=pp*(-gap)*normal[j]*funct1[i];
			}
		}
		// Compute Fc2
		for(int i=0;i<numnode2;i++)
				{
					for(int j=0;j<NDIM;j++)
					{
						fc2[NDIM*i+j]=-pp*(-gap)*normal[j]*funct2[i];
					}
				}
	}
	else	// No penetration --> fc = 0
	{
		// No additional forces --> all entries = 0
		for(int i=0;i<NDIM*numnode1;i++)
			fc1[i]=0;
		for(int i=0;i<NDIM*numnode2;i++)
					fc2[i]=0;
	}
	
//	// Test print of Fc
//	cout<<"************	FD_stiffc ****************"<<endl;
//	for(int i=0;i<NDIM*numnode1;i++)
//	{
//		cout<<"Fc_1_"<<i<<": "<<fc1[i]<<endl;
//	}	
//	for(int i=0;i<NDIM*numnode2;i++)
//	{
//		cout<<"Fc_2_"<<i<<": "<<fc2[i]<<endl;
//	}
	
	// FD-Checks are beginning here
	// loop over ele1pos_
	for(int j=0;j<numnode1;j++)
	{
		for(int i=0;i<NDIM;i++)
		{
			// step forward
			ele1pos_(i,j) += eps; 
//			cout<<"ele1pos_:"<<endl<<ele1pos_<<endl;
						
			// do the closest point projection for the new set
			ClosestPointProjection(FD_XiContact,FD_elementscolinear);
			
			// Call function to fill variables for shape functions and their derivatives
			GetShapeFunctions(FD_funct1,FD_funct2,FD_deriv1,FD_deriv2,FD_secondderiv1,FD_secondderiv2,FD_XiContact);
	
			// Call function to fill variabeles with coordinates and derivatives of the two contacting points
			ComputeCoordsAndDerivs(FD_x1,FD_x2,FD_dx1,FD_dx2,FD_ddx1,FD_ddx2,FD_funct1,FD_funct2,FD_deriv1,FD_deriv2,FD_secondderiv1,FD_secondderiv2,numnode1,numnode2);
	
			// Call function to compute scaled normal and gap in possible contact point
			ComputeNormal(FD_normal,FD_gap,FD_norm,FD_x1,FD_x2);
					
			// local vector for contact forces from FD
			vector<double> FD_fc1(NDIM*numnode1);
			vector<double> FD_fc2(NDIM*numnode2);
						
			// Compute FD_Fc1
			for(int k=0;k<numnode1;k++)
			{
				for(int m=0;m<NDIM;m++)
				{
					FD_fc1[NDIM*k+m]= pp*(-FD_gap)*FD_normal[m]*FD_funct1[k];
					FD_fc2[NDIM*k+m]=-pp*(-FD_gap)*FD_normal[m]*FD_funct2[k];
				}
			}
			
//			// Test print of FD_fc1 and FD_fc2
//			for(int k=0;k<NDIM*numnode1;k++)
//			{
//				cout<<"FD_fc_1_"<<k<<": "<<FD_fc1[k]<<endl;
//			}	
//			for(int k=0;k<NDIM*numnode2;k++)
//			{
//				cout<<"FD_fc_2_"<<k<<": "<<FD_fc2[k]<<endl;
//			}
//			cout<<endl;
			
			// FD-approximation
			for(int k=0;k<NDIM*numnode1;k++)		
			{
				FD_stiffc1(k,NDIM*j+i) = (FD_fc1[k] - fc1[k])/eps;								// left block of FD_stiffc1
				FD_stiffc2(k,NDIM*j+i) = (FD_fc2[k] - fc2[k])/eps;								// left block of FD_stiffc2
			}
										
			// Undo step forward
			ele1pos_(i,j) -= eps;
		}
	}
	// loop over ele2pos_
	for(int j=0;j<numnode2;j++)
	{
		for(int i=0;i<NDIM;i++)
		{
			// step forward
			ele2pos_(i,j) += eps; 
//			cout<<"ele2pos_:"<<endl<<ele2pos_<<endl;
			
			// do the closest point projection for the new set
			ClosestPointProjection(FD_XiContact,FD_elementscolinear);
			
			// Call function to fill variables for shape functions and their derivatives
			GetShapeFunctions(FD_funct1,FD_funct2,FD_deriv1,FD_deriv2,FD_secondderiv1,FD_secondderiv2,FD_XiContact);
	
			// Call function to fill variabeles with coordinates and derivatives of the two contacting points
			ComputeCoordsAndDerivs(FD_x1,FD_x2,FD_dx1,FD_dx2,FD_ddx1,FD_ddx2,FD_funct1,FD_funct2,FD_deriv1,FD_deriv2,FD_secondderiv1,FD_secondderiv2,numnode1,numnode2);
	
			// Call function to compute scaled normal and gap in possible contact point
			ComputeNormal(FD_normal,FD_gap,FD_norm,FD_x1,FD_x2);
			
			// local vector for contact forces from FD
			vector<double> FD_fc1(NDIM*numnode1);
			vector<double> FD_fc2(NDIM*numnode2);
			
			// Compute FD_Fc2
			for(int k=0;k<numnode2;k++)
			{
				for(int m=0;m<NDIM;m++)
				{
//					cout<<"Index Fc_2: "<<NDIM*k+m<<endl;
					FD_fc1[NDIM*k+m]= pp*(-FD_gap)*FD_normal[m]*FD_funct1[k];
					FD_fc2[NDIM*k+m]= -pp*(-FD_gap)*FD_normal[m]*FD_funct2[k];
				}
			}
	
//			// Test print of FD_fc1 and FD_fc2
//			for(int k=0;k<NDIM*numnode1;k++)
//			{
//				cout<<"FD_fc_1_"<<k<<": "<<FD_fc1[k]<<endl;
//			}	
//			for(int k=0;k<NDIM*numnode2;k++)
//			{
//				cout<<"FD_fc_2_"<<k<<": "<<FD_fc2[k]<<endl;
//			}
//			cout<<endl;
			
			// FD-approximation
			for(int k=0;k<NDIM*numnode2;k++)		
			{
				FD_stiffc1(k,NDIM*numnode1+NDIM*j+i) = (FD_fc1[k] - fc1[k])/eps;	// right block of FC_stiffc1
				FD_stiffc2(k,NDIM*numnode1+NDIM*j+i) = (FD_fc2[k] - fc2[k])/eps;	// right block of FD_stiffc2
			}
						
			// Undo step forward
			ele2pos_(i,j) -= eps;
		}
	}
	
	// Print the FD-Approximations to compare with analytical values
//	cout<<endl<<"Finite-Differences-Approximation of stiffc1 and stiffc2:"<<endl<<endl;
//	cout<<"FD_stiffc1:"<<endl<<FD_stiffc1<<endl<<endl<<"FD_stiffc2:"<<endl<<FD_stiffc2<<endl;
	
	// Compare FD-Approximations and analytical values and give a warning, if deviation is to large
	for(int i=0;i<NDIM*numnode1;i++)
	{
		for(int j=0;j<NDIM*(numnode1+numnode2);j++)
		{
			if(stiffc1(i,j) > toltreshold)							// then check relative error
			{
				if(fabs((FD_stiffc1(i,j)-stiffc1(i,j))/stiffc1(i,j)) > rtol)		// compare stiffc1 and FD_stiffc1 elementwise
				{
					cout<<"FD_stiffc1_"<<i<<"_"<<j<<" = "<<FD_stiffc1(i,j)<<"	stiffc1_"<<i<<"_"<<j<<" = ";
					cout<<stiffc1(i,j)<<"	********* WARNING from FD-Check: relative error ";
					cout<<(FD_stiffc1(i,j)-stiffc1(i,j))/stiffc1(i,j)<<" > relative tolerance "<<rtol;
					cout<<" *********"<<endl;
				}
			}
			else																				// check absolute error
			{
				if(fabs(FD_stiffc1(i,j)-stiffc1(i,j)) > atol)		// compare stiffc1 and FD_stiffc1 elementwise
				{
					cout<<"FD_stiffc1_"<<i<<"_"<<j<<" = "<<FD_stiffc1(i,j)<<"	stiffc1_"<<i<<"_"<<j<<" = ";
					cout<<stiffc1(i,j)<<"	********* WARNING from FD-Check: absolute error ";
					cout<<FD_stiffc1(i,j)-stiffc1(i,j)<<" > absolute tolerance "<<atol<<" *********"<<endl;
				}
			}
		}
	}
	for(int i=0;i<NDIM*numnode2;i++)
	{
		for(int j=0;j<NDIM*(numnode1+numnode2);j++)
		{
			if(stiffc2(i,j) > toltreshold)							// then check relative error
			{
				if(fabs((FD_stiffc2(i,j)-stiffc2(i,j))/stiffc2(i,j)) > rtol)		// compare stiffc2 and FD_stiffc2 elementwise
				{
					cout<<"FD_stiffc2_"<<i<<"_"<<j<<" = "<<FD_stiffc2(i,j)<<"	stiffc2_"<<i<<"_"<<j<<" = ";
					cout<<stiffc2(i,j)<<"	********* WARNING from FD-Check: relative error ";
					cout<<(FD_stiffc2(i,j)-stiffc2(i,j))/stiffc2(i,j)<<" > relative tolerance "<<rtol;
					cout<<" *********"<<endl;
				}
			}
			else																				// check absolute error
			{
				if(fabs(FD_stiffc2(i,j)-stiffc2(i,j)) > atol)		// compare stiffc2 and FD_stiffc2 elementwise
				{
					cout<<"FD_stiffc2_"<<i<<"_"<<j<<" = "<<FD_stiffc2(i,j)<<"	stiffc2_"<<i<<"_"<<j<<" = ";
					cout<<stiffc2(i,j)<<"	********* WARNING from FD-Check: absolute error ";
					cout<<FD_stiffc2(i,j)-stiffc2(i,j)<<" > absolute tolerance "<<atol<<" *********"<<endl;
				}
			}
		}
	}
	
	return;
}

#endif  // #ifdef CCADISCRET 
