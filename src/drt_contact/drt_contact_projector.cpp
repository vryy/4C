/*!----------------------------------------------------------------------
\file drt_contact_projector.cpp
\brief A class to perform projections of nodes onto opposing CElements

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_contact_projector.H"
#include "drt_contact_interface.H"
#include "drt_celement.H"
#include "drt_cnode.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 01/08|
 *----------------------------------------------------------------------*/
CONTACT::Projector::Projector(bool twoD) :
twoD_(twoD)
{
}

/*----------------------------------------------------------------------*
 |  Project a node along its nodal normal (public)            popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::Project_NodalNormal(CONTACT::CNode& node,
		 																		     CONTACT::CElement& ele,
		 																		     double xi[])
{
	bool ok = true;
	if (IsTwoDimensional())
	{
#ifdef DEBUG
		// two-dimensional version of the problem
		//cout << "CONTACT::Projector::Project_NodalNormal" << endl
		//		 << "Ready for projection of slave CNode " << node.Id()
		//		 << " onto master CElement " << ele.Id() << endl;
#endif // #ifdef DEBUG
		
		// local Newton iteration for xi, start in the element middle
		double eta[2] = {0.0, 0.0};
		double deta = 0.0;
		double F = 0.0;
		double dF = 0.0;
		int k=0;
		
		for (k=0;k<CONTACT_MAXITER;++k)
		{
			F=Evaluate_F_NodalNormal(node,ele,eta);
			if (abs(F) < CONTACT_CONVTOL) break;
			dF=Evaluate_gradF_NodalNormal(node,ele,eta);
			deta=(-F)/dF;
			eta[0]+=deta;
		}
		
		// get the result
		xi[0]=eta[0];
		
		// Newton iteration unconverged
		if (abs(F) > CONTACT_CONVTOL)
		{
			ok = false;
			xi[0] = 9999.99;
			
			// Here (S->M projection) we only give a warning, no error!!!
			// This iteration sometimes diverges, when the projection is far off.
			// Although these cases are harmless, as these nodes then do not participate in
			// the overlap detection anyway one should check these warnings (at the moment by hand)!
			// (FIXME: Automatic check, if unconverged projections are really harmless!)
			cout << "***WARNING*** Project_NodalNormal:" << endl << "Newton unconverged for NodeID "
					 << node.Id() << " and CElementID " << ele.Id() << endl;	
		}
		
#ifdef DEBUG			
		// Newton iteration converged
		else
		{
			cout << "Newton iteration converged in " << k << " step(s)!" << endl
			  	 << "The result is: " << xi[0] << endl;
		}
#endif // #ifdef DEBUG
		
	} // if (IsTwoDimesional())
	
	else
	{
		// three-dimensional version of the problem
		ok = false;
		dserror("ERROR: Project_NodalNormal: 3D version not yet implemented!");
	}
	
	return true;
}

/*----------------------------------------------------------------------*
 |  Project a node along element's normal field (public)      popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::Project_ElementNormal(CONTACT::CNode& node,
		 																			     CONTACT::CElement& ele,
		 																			     double xi[])
{
	bool ok = true;
	
	if (IsTwoDimensional())
	{
#ifdef DEBUG
		// two-dimensional version of the problem
		//cout << "CONTACT::Projector::Project_ElementNormal" << endl
		//		 << "Ready for projection of master CNode " << node.Id()
		//		 << " onto slave CElement " << ele.Id() << endl;
#endif // #ifdef DEBUG
			
		// local Newton iteration for xi, start in the element middle
		double eta[2] = {0.0, 0.0};
		double deta = 0.0;
		double F = 0.0;
		double dF = 0.0;
		int k=0;
		
		for (k=0;k<CONTACT_MAXITER;++k)
		{
			F=Evaluate_F_ElementNormal(node,ele,eta);
			if (abs(F) < CONTACT_CONVTOL) break;
			dF=Evaluate_gradF_ElementNormal(node,ele,eta);
			deta=(-F)/dF;
			eta[0]+=deta;
		}
		
		// get the result
		xi[0]=eta[0];
		
		// Newton iteration unconverged
		if (abs(F) > CONTACT_CONVTOL)
		{
			ok = false;
			xi[0] = 9999.99;
			
			// Here (M->S projection) we only give a warning, no error!!!
			// This iteration sometimes diverges, when the projection is far off.
			// Although these cases are harmless, as these nodes then do not participate in
			// the overlap detection anyway one should check these warnings (at the moment by hand)!
			// (FIXME: Automatic check, if unconverged projections are really harmless!)
			cout << "***WARNING*** Project_ElementNormal:" << endl << "Newton unconverged for NodeID "
					 << node.Id() << " and CElementID " << ele.Id() << endl;
		}
		
#ifdef DEBUG			
		// Newton iteration converged
		else
		{
			cout << "Newton iteration converged in " << k << " step(s)!" << endl
			  	 << "The result is: " << xi[0] << endl;
		}
#endif // #ifdef DEBUG
		
	} // if (IsTwoDimesional())
		
	else
	{
		// three-dimensional version of the problem
		ok = false;
		dserror("ERROR: Project_ElementNormal: 3D version not yet implemented!");
	}
		
	return ok;
}

/*----------------------------------------------------------------------*
 |  Project a Gauss point along its normal (public)           popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::Projector::Project_GaussPoint(CONTACT::CElement& gpele,
																						const double* gpeta,
		 																		    CONTACT::CElement& ele,
		 																		    double xi[])
{
	bool ok = true;
	if (IsTwoDimensional())
	{
		// two-dimensional version of the problem
#ifdef DEBUG
		//cout << "CONTACT::Projector::Project_GaussPoint" << endl
		//		 << "Ready for projection of GP at xi=" << gpeta[0]
		//		 << " from slave Celement " << gpele.Id()
		//		 << " onto master CElement " << ele.Id() << endl;
#endif // #ifdef DEBUG
		
		// collect necessary data (slave side, for GP)
		int nnodes = gpele.NumNode();
		vector<double> val(nnodes);
		vector<double> deriv(nnodes);
		LINALG::SerialDenseMatrix coord(3,nnodes);
		DRT::Node** mynodes = gpele.Nodes();
		if(!mynodes) dserror("ERROR: Project_GaussPoint: Null pointer!");
				
		// get shape function values and derivatives at gpeta
		ele.EvaluateShape_1D(gpeta, val, deriv, nnodes);

		// get interpolated GP normal and GP coordinates
		double gpn[3] = {0.0, 0.0, 0.0};
		double gpx[3] = {0.0, 0.0, 0.0};
		for (int i=0;i<nnodes;++i)
		{
			CNode* mycnode = static_cast<CNode*> (mynodes[i]);
			
			gpn[0]+=val[i]*mycnode->n()[0];
			gpn[1]+=val[i]*mycnode->n()[1];
			gpn[2]+=val[i]*mycnode->n()[2];
			
			coord(0,i) = mycnode->xspatial()[0];
			coord(1,i) = mycnode->xspatial()[1];
			coord(2,i) = mycnode->xspatial()[2];
			
			gpx[0]+=val[i]*coord(0,i);
			gpx[1]+=val[i]*coord(1,i);
			gpx[2]+=val[i]*coord(2,i);
		}
		
		// local Newton iteration for xi, start in the element middle
		double eta[2] = {0.0, 0.0};
		double deta = 0.0;
		double F = 0.0;
		double dF = 0.0;
		int k=0;
		
		for (k=0;k<CONTACT_MAXITER;++k)
		{
			F=Evaluate_F_GaussPoint(gpx,gpn,ele,eta);
			if (abs(F) < CONTACT_CONVTOL) break;
			dF=Evaluate_gradF_GaussPoint(gpn,ele,eta);
			deta=(-F)/dF;
			eta[0]+=deta;
		}
		
		// Newton iteration unconverged
		if (abs(F) > CONTACT_CONVTOL)
		{
			ok = false;
			dserror("ERROR: Project_GaussPoint: Newton unconverged for GP at xi=%"
							" from CElementID %", gpeta[0], gpele.Id());
		}
		
		// Newton iteration converged
		xi[0]=eta[0];
#ifdef DEBUG
		cout << "Newton iteration converged in " << k << " step(s)!" << endl
				 << "The result is: " << xi[0] << endl;
#endif // #ifdef DEBUG
	}
	
	else
	{
		// three-dimensional version of the problem
		ok = false;
		dserror("ERROR: Project_GaussPoint: 3D version not yet implemented!");
	}
	
	return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for nodal normal case (public)                 popp 01/08|
 *----------------------------------------------------------------------*/
double CONTACT::Projector::Evaluate_F_NodalNormal(CONTACT::CNode& node,
     																							CONTACT::CElement& ele,
     																							const double* eta)
{
	/* Evaluate the function F(eta) = ( Ni * xim - xs ) x ns,
		 or to be more precise the third component of this vector function!
		
		   Ni  shape functions of element
		   xim coords of element nodes (master side)
		   xs  coords of node to be projected (slave side)
		   ns	 outward normal of node to be projected (slave side)          */
	
	// collect necessary data (master side)
	int nnodes = ele.NumNode();
	vector<double> val(nnodes);
	vector<double> deriv(nnodes);
	
	// get shape function values and derivatives at eta
	ele.EvaluateShape_1D(eta, val, deriv, nnodes);

	// build interpolation of master node coordinates for current eta
	double Nx[3] = {0.0, 0.0, 0.0};
	ele.LocalToGlobal(eta,Nx,true);
	
	// subtract slave node coordinates
	Nx[0]-=node.xspatial()[0];
	Nx[1]-=node.xspatial()[1];
	Nx[2]-=node.xspatial()[2];

	// calculate F
	return (Nx[0]*node.n()[1]-Nx[1]*node.n()[0]);
}

/*----------------------------------------------------------------------*
 |  Evaluate gradF for nodal normal case (public)             popp 01/08|
 *----------------------------------------------------------------------*/
double CONTACT::Projector::Evaluate_gradF_NodalNormal(CONTACT::CNode& node,
     																							    CONTACT::CElement& ele,
     																							    const double* eta)
{
	/* Evaluate the function gradF(eta)
	   = Ni,eta * xim * nys - Ni,eta * yim * nxs,
		
		   Ni,eta    shape function derivatives of element
		   xim, yim  coords of element nodes (master side)
		   nxs, nys	 outward normal of node to be projected (slave side)   */
	
	// collect necessary data (master side)
	int nnodes = ele.NumNode();
	vector<double> val(nnodes);
	vector<double> deriv(nnodes);
			
	// get shape function values and derivatives at eta
	ele.EvaluateShape_1D(eta, val, deriv, nnodes);

	// build interpolation of master node coordinates for current eta
	// use shape function derivatives for interpolation
	double Nxeta[3] = {0.0, 0.0, 0.0};
	ele.LocalToGlobal(eta,Nxeta,false);
	
	// calculate gradF
	return (Nxeta[0]*node.n()[1]-Nxeta[1]*node.n()[0]);
}

/*----------------------------------------------------------------------*
 |  Evaluate F for element normal case (public)               popp 01/08|
 *----------------------------------------------------------------------*/
double CONTACT::Projector::Evaluate_F_ElementNormal(CONTACT::CNode& node,
     																							  CONTACT::CElement& ele,
     																							  const double* eta)
{
	/* Evaluate the function F(eta) = ( Ni * xis - xm ) x ( Nj * njs),
	   or to be more precise the third component of this vector function!
			
		   Ni  shape functions of element
		   xis coords of element nodes (slave side)
		   xm  coords of node to be projected (master side)
		   nis outward normals of element nodes (slave side)                */
		
	// collect necessary data (slave side)
	int nnodes = ele.NumNode();
	vector<double> val(nnodes);
	vector<double> deriv(nnodes);
	LINALG::SerialDenseMatrix coord(3,nnodes);
	DRT::Node** mynodes = ele.Nodes();
	if(!mynodes) dserror("ERROR: Evaluate_F_ElementNormal: Null pointer!");
			
	// get shape function values and derivatives at eta
	ele.EvaluateShape_1D(eta, val, deriv, nnodes);

	// get interpolated normal and proj. coordinates for current eta
	double Nn[3] = {0.0, 0.0, 0.0};
	double Nx[3] = {0.0, 0.0, 0.0};
	for (int i=0;i<nnodes;++i)
	{
		CNode* mycnode = static_cast<CNode*> (mynodes[i]);
		Nn[0]+=val[i]*mycnode->n()[0];
		Nn[1]+=val[i]*mycnode->n()[1];
		Nn[2]+=val[i]*mycnode->n()[2];
		
		coord(0,i) = mycnode->xspatial()[0];
		coord(1,i) = mycnode->xspatial()[1];
		coord(2,i) = mycnode->xspatial()[2];
		
		Nx[0]+=val[i]*coord(0,i);
		Nx[1]+=val[i]*coord(1,i);
		Nx[2]+=val[i]*coord(2,i);
	}
		
	// subtract master node coordinates
	Nx[0]-=node.xspatial()[0];
	Nx[1]-=node.xspatial()[1];
	Nx[2]-=node.xspatial()[2];
	
	// calculate F
	return (Nx[0]*Nn[1]-Nx[1]*Nn[0]);
}

/*----------------------------------------------------------------------*
 |  Evaluate gradF for element normal case (public)           popp 01/08|
 *----------------------------------------------------------------------*/
double CONTACT::Projector::Evaluate_gradF_ElementNormal(CONTACT::CNode& node,
     																							      CONTACT::CElement& ele,
     																							      const double* eta)
{
	/* Evaluate the function gradF(eta)
	 	 = ( Ni,eta * xis ) * ( Nj * nyjs )
	 	 + ( Ni * xis - xm ) * ( Nj,eta * nyjs )
	 	 - ( Ni,eta * yis ) * ( Nj * nxjs )
	 	 - ( Ni * yis - ym ) * ( Nj,eta * nxjs )
				
		   Ni,eta  		shape function derivatives of element
		   xis, yis   coords of element nodes (slave side)
		   xm, ym     coords of node to be projected (master side)
		   nxjs, nyjs outward normals of element nodes (slave side)         */
			
	// collect necessary data (slave side)
	int nnodes = ele.NumNode();
	vector<double> val(nnodes);
	vector<double> deriv(nnodes);
	LINALG::SerialDenseMatrix coord(3,nnodes);
	DRT::Node** mynodes = ele.Nodes();
	if(!mynodes) dserror("ERROR: Evaluate_gradF_ElementNormal: Null pointer!");
				
	// get shape function values and derivatives at eta
	ele.EvaluateShape_1D(eta, val, deriv, nnodes);

	// get interpolated normal and proj. coordinates for current eta
	double Nn[3] = {0.0, 0.0, 0.0};
	double Nneta[3] = {0.0, 0.0, 0.0};
	double Nx[3] = {0.0, 0.0, 0.0};
	double Nxeta[3] = {0.0, 0.0, 0.0};
	for (int i=0;i<nnodes;++i)
	{
		CNode* mycnode = static_cast<CNode*> (mynodes[i]);
		
		Nn[0]+=val[i]*mycnode->n()[0];
		Nn[1]+=val[i]*mycnode->n()[1];
		Nn[2]+=val[i]*mycnode->n()[2];
		
		Nneta[0]+=deriv[i]*mycnode->n()[0];
		Nneta[1]+=deriv[i]*mycnode->n()[1];
		Nneta[2]+=deriv[i]*mycnode->n()[2];
			
		coord(0,i) = mycnode->xspatial()[0];
		coord(1,i) = mycnode->xspatial()[1];
		coord(2,i) = mycnode->xspatial()[2];
		
		Nx[0]+=val[i]*coord(0,i);
		Nx[1]+=val[i]*coord(1,i);
		Nx[2]+=val[i]*coord(2,i);
				
		Nxeta[0]+=deriv[i]*coord(0,i);
		Nxeta[1]+=deriv[i]*coord(1,i);
		Nxeta[2]+=deriv[i]*coord(2,i);
	}
			
	// subtract master node coordinates
	Nx[0]-=node.xspatial()[0];
	Nx[1]-=node.xspatial()[1];
	Nx[2]-=node.xspatial()[2];
		
	// calculate gradF
	double gradF =   Nxeta[0]*Nn[1] + Nx[0]*Nneta[1]
	               - Nxeta[1]*Nn[0] - Nx[1]*Nneta[0];
	return gradF;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for Gauss point case (public)                  popp 01/08|
 *----------------------------------------------------------------------*/
double CONTACT::Projector::Evaluate_F_GaussPoint(const double* gpx,
																								 const double* gpn,
     																						 CONTACT::CElement& ele,
     																						 const double* eta)
{
	/* Evaluate the function F(eta) = ( Ni * xim - gpx ) x gpn,
		 or to be more precise the third component of this vector function!
		
		   Ni  shape functions of element (master side)
		   xim coords of element nodes (master side)
		   gpx coords of GP to be projected (slave side)
		   gpn outward normal of GP to be projected (slave side)          */
	
	// collect necessary data (master side)
	int nnodes = ele.NumNode();
	vector<double> val(nnodes);
	vector<double> deriv(nnodes);
		
	// get shape function values and derivatives at eta
	ele.EvaluateShape_1D(eta, val, deriv, nnodes);

	// build interpolation of master node coordinates for current eta
	double Nx[3] = {0.0, 0.0, 0.0};
	ele.LocalToGlobal(eta,Nx,true);
	
	// subtract GP coordinates
	Nx[0]-=gpx[0];
	Nx[1]-=gpx[1];
	Nx[2]-=gpx[2];

	// calculate F
	return (Nx[0]*gpn[1]-Nx[1]*gpn[0]);
}

/*----------------------------------------------------------------------*
 |  Evaluate gradF for Gauss point case (public)              popp 01/08|
 *----------------------------------------------------------------------*/
double CONTACT::Projector::Evaluate_gradF_GaussPoint(const double* gpn,
																								 	   CONTACT::CElement& ele,
																								 	   const double* eta)
{
	/* Evaluate the function gradF(eta)
	 	 = Ni,eta * xim * gpny - Ni,eta * yim * gpnx,
		
		   Ni,eta     shape function derivatives of element (master side)
		   xim, yim   coords of element nodes (master side)
		   gpnx, gpny outward normal of GP to be projected (slave side)   */
	
	// collect necessary data (master side)
	int nnodes = ele.NumNode();
	vector<double> val(nnodes);
	vector<double> deriv(nnodes);
			
	// get shape function values and derivatives at eta
	ele.EvaluateShape_1D(eta, val, deriv, nnodes);

	// build interpolation of master node coordinates for current eta
	// use shape fuvntion derivatives for interpolation
	double Nxeta[3] = {0.0, 0.0, 0.0};
	ele.LocalToGlobal(eta,Nxeta,false);

	// calculate gradF
	return (Nxeta[0]*gpn[1]-Nxeta[1]*gpn[0]);
}

#endif //#ifdef CCADISCRET
