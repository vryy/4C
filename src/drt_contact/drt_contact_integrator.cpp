/*!----------------------------------------------------------------------
\file drt_contact_integrator.cpp
\brief A class to perform intgerations of Mortar matrices on the overlap
 			 of two CElements in 1D and 2D

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_contact_integrator.H"
#include "drt_celement.H"
#include "drt_contact_projector.H"
#include "../drt_lib/drt_utils.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 01/08|
 *----------------------------------------------------------------------*/
CONTACT::Integrator::Integrator(int ngp, bool oneD) :
oneD_(oneD),
ngp_(ngp)
{
	if (oneD)
	{
		coords_.resize(ngp_);
		weights_.resize(ngp_);
		
		switch(ngp_)
		{
		case 1:
		{
			const DRT::Utils::IntegrationPoints1D intpoints(DRT::Utils::intrule_line_1point);
			for (int i=0;i<intpoints.nquad;++i)
			{
				coords_[i]=intpoints.qxg[i];
				weights_[i]=intpoints.qwgt[i];
			}
			break;
		}
		case 2:
		{
			const DRT::Utils::IntegrationPoints1D intpoints(DRT::Utils::intrule_line_2point);
			for (int i=0;i<intpoints.nquad;++i)
			{
				coords_[i]=intpoints.qxg[i];
				weights_[i]=intpoints.qwgt[i];
			}
			break;
		}
		case 3:
		{
			const DRT::Utils::IntegrationPoints1D intpoints(DRT::Utils::intrule_line_3point);
			for (int i=0;i<intpoints.nquad;++i)
			{
				coords_[i]=intpoints.qxg[i];
				weights_[i]=intpoints.qwgt[i];
			}
			break;
		}
		case 4:
		{
			const DRT::Utils::IntegrationPoints1D intpoints(DRT::Utils::intrule_line_4point);
			for (int i=0;i<intpoints.nquad;++i)
			{
				coords_[i]=intpoints.qxg[i];
				weights_[i]=intpoints.qwgt[i];
			}
			break;
		}
		case 5:
		{
			const DRT::Utils::IntegrationPoints1D intpoints(DRT::Utils::intrule_line_5point);
			for (int i=0;i<intpoints.nquad;++i)
			{
				coords_[i]=intpoints.qxg[i];
				weights_[i]=intpoints.qwgt[i];
			}
			break;
		}	
		default:
			dserror("ERROR: Integrator: Given no. of Gauss points not implemented!");
		} // switch (ngp_)
		
		
	} // if (oneD)
	
	else
		dserror("ERROR: Integrator: 2D case not yet implemented!");
}

/*----------------------------------------------------------------------*
 |  Integrate a 1D slave element overlap                      popp 01/08|
 |	This method integrates 2 functions on the same (slave) CEelement    |
 |	from given local coordinates sxia to sxib														|
 |	Output is an Epetra_SerialDenseMatrix holding the int. values 			|
 *----------------------------------------------------------------------*/
RCP<Epetra_SerialDenseMatrix> CONTACT::Integrator::Integrate_D(CONTACT::CElement& sele,
																					 										double sxia, double sxib)
{
	//check input data
	if (!sele.IsSlave())
		dserror("ERROR: Integrate_D called on a non-slave CElement!");
	if ((sxia<-1.0) || (sxib>1.0))
		dserror("ERROR: Integrate_D called with infeasible slave limits!");
	
	// create empty Ddense object and wrap it with RCP
	int nrow = sele.NumNode();
	int ncol = nrow;
	RCP<Epetra_SerialDenseMatrix> Ddense = rcp(new Epetra_SerialDenseMatrix(nrow,ncol));
	
	// create empty vectors for shape fct. evaluation
	vector<double> val(nrow);
	vector<double> deriv(nrow);
	vector<double> dualval(nrow);
	vector<double> dualderiv(nrow);
	
	// get nodal coords for Jacobian evaluation
	LINALG::SerialDenseMatrix coord(3,nrow);
	coord = sele.GetNodalCoords();
	
	// loop over all Gauss points for integration
	for (int gp=0;gp<nGP();++gp)
	{
		double eta[2] = {0.0, 0.0};
		eta[0] = Coordinate(gp);
		double wgt = Weight(gp);
		
		// coordinate transformation sxi->eta (CElement->Overlap)
		double sxi[2] = {0.0, 0.0};
		sxi[0] = 0.5*(1-eta[0])*sxia + 0.5*(1+eta[0])*sxib;
		
		// evaluate trace space and dual space shape functions
		sele.EvaluateShape_1D(sxi,val,deriv,nrow);
		sele.EvaluateShape_Dual1D(sxi,dualval,dualderiv,nrow);
		
		// evaluate the two Jacobians
		double dxdsxi = sele.Jacobian_1D(val,deriv,coord);
		double dsxideta = -0.5*sxia + 0.5*sxib;
		
		/* loop over all Ddense matrix entries
		   nrow represents the dofs !!!
		   ncol represents the Lagrange multipliers !!!
		   (although this does not really matter here for Ddense,
		   as it will turn out to be diagonal anyway)             */
		for (int j=0;j<nrow;++j)
		{
			for (int k=0;k<ncol;++k)
			{
				// multiply the two shape functions
				double prod = val[j]*dualval[k];
				// add current Gauss point's contribution to Ddense  
				(*Ddense)(j,k) += prod*dxdsxi*dsxideta*wgt; 
			}
		}	
	} // for (int gp=0;gp<nGP();++gp)
	
#ifdef DEBUG
	cout << endl << *Ddense;
#endif // #ifdef DEBUG
	
	return Ddense;
}

/*----------------------------------------------------------------------*
 |  Integrate a 1D slave / master overlap                     popp 01/08|
 |	This method integrates a slave side function (dual shape fct.)      |
 |  and a master side function (standard shape fct.) from given local   |
 |	coordinates sxia to sxib. The corresponding master side local       |
 |  element coordinates given by mxia and mxib												  |
 |	Output is an Epetra_SerialDenseMatrix holding the int. values 			|
 *----------------------------------------------------------------------*/
RCP<Epetra_SerialDenseMatrix> CONTACT::Integrator::Integrate_M(CONTACT::CElement& sele,
																					 										double sxia, double sxib,
																					 										CONTACT::CElement& mele,
																					 										double mxia, double mxib)
{
	// check input data
	if ((!sele.IsSlave()) || (mele.IsSlave()))
		dserror("ERROR: Integrat_eM called on a wrong type of CElement pair!");
	if ((sxia<-1.0) || (sxib>1.0))
		dserror("ERROR: Integrate_M called with infeasible slave limits!");
	if ((mxia<-1.0) || (mxib>1.0))
			dserror("ERROR: Integrate_M called with infeasible master limits!");
	
	// create empty Mdense object and wrap it with RCP
	int nrow = sele.NumNode();
	int ncol = mele.NumNode();
	RCP<Epetra_SerialDenseMatrix> Mdense = rcp(new Epetra_SerialDenseMatrix(nrow,ncol));
	
	// create empty vectors for shape fct. evaluation
	vector<double> sval(nrow);
	vector<double> sderiv(nrow);
	vector<double> mval(ncol);
	vector<double> mderiv(ncol);
	vector<double> dualval(nrow);
	vector<double> dualderiv(nrow);

	// get slave nodal coords for Jacobian evaluation
	LINALG::SerialDenseMatrix scoord(3,nrow);
	scoord = sele.GetNodalCoords();
	
	// loop over all Gauss points for integration
	for (int gp=0;gp<nGP();++gp)
	{
		double eta[2] = {0.0, 0.0};
		eta[0] = Coordinate(gp);
		double wgt = Weight(gp);
		
		
		double sxi[2] = {0.0, 0.0};
		double mxi[2] = {0.0, 0.0};
		
		// coordinate transformation sxi->eta (slave CElement->Overlap)
		sxi[0] = 0.5*(1-eta[0])*sxia + 0.5*(1+eta[0])*sxib;
		
		// project Gauss point onto master element
		CONTACT::Projector projector(true);
		projector.Project_GaussPoint(sele,sxi,mele,mxi);
		
		// check GP projection
		if ((mxi[0]<mxia) || (mxi[0]>mxib))
			dserror("ERROR: Integrate_M: Gauss point projection failed!");
		
		// evaluate dual space shape functions (on slave element)
		sele.EvaluateShape_Dual1D(sxi,dualval,dualderiv,nrow);
		
		// evaluate trace space shape functions (on both elements)
		sele.EvaluateShape_1D(sxi,sval,sderiv,nrow);
		mele.EvaluateShape_1D(mxi,mval,mderiv,ncol);
		
		// evaluate the two slave side Jacobians
		double dxdsxi = sele.Jacobian_1D(sval,sderiv,scoord);
		double dsxideta = -0.5*sxia + 0.5*sxib;
		
		/* loop over all Mdense matrix entries
		   nrow represents the slave Lagrange multipliers !!!
		   ncol represents the master dofs !!!
		   (this DOES matter here for Mdense, as it might
		   sometimes be rectangular, not quadratic!)              */
		for (int j=0;j<nrow;++j)
		{
			for (int k=0;k<ncol;++k)
			{
				// multiply the two shape functions
				double prod = dualval[j]*mval[k];
				// add current Gauss point's contribution to Mdense  
				(*Mdense)(j,k) += prod*dxdsxi*dsxideta*wgt; 
			}
		}	
	} // for (int gp=0;gp<nGP();++gp)
	
#ifdef DEBUG
	cout << endl << *Mdense;
#endif // #ifdef DEBUG
	
	return Mdense;
}

/*----------------------------------------------------------------------*
 |  Integrate gap on a 1D slave / master overlap              popp 01/08|
 |	This method integrates a slave side function (dual shape fct.)      |
 |  and the gap function g = ( ( sx - mx ) * n )  from given local      |
 |	coordinates sxia to sxib. The corresponding master side local       |
 |  element coordinates given by mxia and mxib												  |
 |	Output is an Epetra_SerialDenseVector holding the int. values 			|
 *----------------------------------------------------------------------*/
RCP<Epetra_SerialDenseVector> CONTACT::Integrator::Integrate_g(CONTACT::CElement& sele,
																					 										double sxia, double sxib,
																					 										CONTACT::CElement& mele,
																					 										double mxia, double mxib)
{
	// check input data
	if ((!sele.IsSlave()) || (mele.IsSlave()))
		dserror("ERROR: Integrate_g called on a wrong type of CElement pair!");
	if ((sxia<-1.0) || (sxib>1.0))
		dserror("ERROR: Integrate_g called with infeasible slave limits!");
	if ((mxia<-1.0) || (mxib>1.0))
			dserror("ERROR: Integrate_g called with infeasible master limits!");
	
	// create empty gdense object and wrap it with RCP
	int nrow = sele.NumNode();
	int ncol = mele.NumNode();
	RCP<Epetra_SerialDenseVector> gdense = rcp(new Epetra_SerialDenseVector(nrow));
	
	// create empty vectors for shape fct. evaluation
	vector<double> sval(nrow);
	vector<double> sderiv(nrow);
	vector<double> mval(ncol);
	vector<double> mderiv(ncol);
	vector<double> dualval(nrow);
	vector<double> dualderiv(nrow);

	// get slave and master nodal coords for Jacobian / GP evaluation
	LINALG::SerialDenseMatrix scoord(3,nrow);
	LINALG::SerialDenseMatrix mcoord(3,ncol);
	scoord = sele.GetNodalCoords();
	mcoord = mele.GetNodalCoords();
	
	// get slave element nodes themselves for normal evaluation
	DRT::Node** mynodes = sele.Nodes();
	if(!mynodes) dserror("ERROR: Integrate_g: Null pointer!");
	
	// loop over all Gauss points for integration
	for (int gp=0;gp<nGP();++gp)
	{
		double eta[2] = {0.0, 0.0};
		eta[0] = Coordinate(gp);
		double wgt = Weight(gp);
		
		double sxi[2] = {0.0, 0.0};
		double mxi[2] = {0.0, 0.0};
		
		// coordinate transformation sxi->eta (slave CElement->Overlap)
		sxi[0] = 0.5*(1-eta[0])*sxia + 0.5*(1+eta[0])*sxib;
		
		// project Gauss point onto master element
		CONTACT::Projector projector(true);
		projector.Project_GaussPoint(sele,sxi,mele,mxi);
		
		// check GP projection
		if ((mxi[0]<mxia) || (mxi[0]>mxib))
			dserror("ERROR: Integrate_g: Gauss point projection failed!");
		
		// evaluate dual space shape functions (on slave element)
		sele.EvaluateShape_Dual1D(sxi,dualval,dualderiv,nrow);
		
		// evaluate trace space shape functions (on both elements)
		sele.EvaluateShape_1D(sxi,sval,sderiv,nrow);
		mele.EvaluateShape_1D(mxi,mval,mderiv,ncol);
		
		// build interpolation of slave GP normal and coordinates
		double gpn[3] = {0.0,0.0,0.0};
		double sgpx[3] = {0.0, 0.0, 0.0};
		for (int i=0;i<nrow;++i)
		{
			CNode* mycnode = static_cast<CNode*> (mynodes[i]);
			gpn[0]+=sval[i]*mycnode->n()[0];
			gpn[1]+=sval[i]*mycnode->n()[1];
			gpn[2]+=sval[i]*mycnode->n()[2];
						
			sgpx[0]+=sval[i]*scoord(0,i);
			sgpx[1]+=sval[i]*scoord(1,i);
			sgpx[2]+=sval[i]*scoord(2,i);
		}
		
		// normalize interpolated GP normal back to length 1.0 !!!
		double length = sqrt(gpn[0]*gpn[0]+gpn[1]*gpn[1]+gpn[2]*gpn[2]);
		for (int i=0;i<3;++i)
			gpn[i]/=length;
		
		// build interpolation of slave GP coordinates
		double mgpx[3] = {0.0, 0.0, 0.0};
		for (int i=0;i<ncol;++i)
		{
			mgpx[0]+=mval[i]*mcoord(0,i);
			mgpx[1]+=mval[i]*mcoord(1,i);
			mgpx[2]+=mval[i]*mcoord(2,i);
		}
		
		// build normal gap at current GP
		double gap = 0.0;
		for (int i=0;i<3;++i)
			gap+=(mgpx[i]-sgpx[i])*gpn[i];
				
#ifdef DEBUG
		gap = sqrt(gap);
		cout << "GP gap: " << gap << endl;
#endif // #ifdef DEBUG
		
		// evaluate the two slave side Jacobians
		double dxdsxi = sele.Jacobian_1D(sval,sderiv,scoord);
		double dsxideta = -0.5*sxia + 0.5*sxib;
		
		/* loop over all gdense vector entries
			 nrow represents the slave side dofs !!!                */
		for (int j=0;j<nrow;++j)
		{
			double prod = dualval[j]*gap;
			// add current Gauss point's contribution to gdense  
			(*gdense)(j) += prod*dxdsxi*dsxideta*wgt; 
		}
		
	} // for (int gp=0;gp<nGP();++gp)
	
#ifdef DEBUG
	cout << endl << *gdense;
#endif // #ifdef DEBUG
	
	return gdense;
}


#endif //#ifdef CCADISCRET
