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
			const DRT::UTILS::IntegrationPoints1D intpoints(DRT::UTILS::intrule_line_1point);
			for (int i=0;i<intpoints.nquad;++i)
			{
				coords_[i]=intpoints.qxg[i];
				weights_[i]=intpoints.qwgt[i];
			}
			break;
		}
		case 2:
		{
			const DRT::UTILS::IntegrationPoints1D intpoints(DRT::UTILS::intrule_line_2point);
			for (int i=0;i<intpoints.nquad;++i)
			{
				coords_[i]=intpoints.qxg[i];
				weights_[i]=intpoints.qwgt[i];
			}
			break;
		}
		case 3:
		{
			const DRT::UTILS::IntegrationPoints1D intpoints(DRT::UTILS::intrule_line_3point);
			for (int i=0;i<intpoints.nquad;++i)
			{
				coords_[i]=intpoints.qxg[i];
				weights_[i]=intpoints.qwgt[i];
			}
			break;
		}
		case 4:
		{
			const DRT::UTILS::IntegrationPoints1D intpoints(DRT::UTILS::intrule_line_4point);
			for (int i=0;i<intpoints.nquad;++i)
			{
				coords_[i]=intpoints.qxg[i];
				weights_[i]=intpoints.qwgt[i];
			}
			break;
		}
		case 5:
		{
			const DRT::UTILS::IntegrationPoints1D intpoints(DRT::UTILS::intrule_line_5point);
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
	
	// create empty D_seg object and wrap it with RCP
	int nrow = sele.NumNode();
	int ndof = 2;              // up to now we only consider 2D problems!!!
	int ncol = nrow;
	
	RCP<Epetra_SerialDenseMatrix> D_temp = rcp(new Epetra_SerialDenseMatrix(nrow,ncol));
	RCP<Epetra_SerialDenseMatrix> D_seg = rcp(new Epetra_SerialDenseMatrix(nrow*ndof,ncol*ndof));
	
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
		
		/* loop over all D_temp matrix entries
		   nrow represents the dofs !!!
		   ncol represents the Lagrange multipliers !!!
		   (although this does not really matter here for D_seg,
		   as it will turn out to be diagonal anyway)             */
		for (int j=0;j<nrow;++j)
		{
			for (int k=0;k<ncol;++k)
			{
				// multiply the two shape functions
				double prod = val[j]*dualval[k];
				// add current Gauss point's contribution to D_temp  
				(*D_temp)(j,k) += prod*dxdsxi*dsxideta*wgt; 
			}
		}	
	} // for (int gp=0;gp<nGP();++gp)

	// fill D_seg matrix with D_temp matrix entries
	// (each D_temp value is multiplied with a (dof)-unit-matrix)
	for (int j=0;j<nrow*ndof;++j)
	{
		for (int k=0;k<ncol*ndof;++k)
		{
			int jindex = (int)(j/ndof);
			int kindex = (int)(k/ndof);
			// isolate the D_seg entries to be filled
			if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
				(*D_seg)(j,k) = (*D_temp)(jindex,kindex);
		}
	}

	return D_seg;
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
		dserror("ERROR: Integrate_M called on a wrong type of CElement pair!");
	if ((sxia<-1.0) || (sxib>1.0))
		dserror("ERROR: Integrate_M called with infeasible slave limits!");
	if ((mxia<-1.0) || (mxib>1.0))
			dserror("ERROR: Integrate_M called with infeasible master limits!");
	
	// create empty M_seg object and wrap it with RCP
	int nrow = sele.NumNode();
	int ncol = mele.NumNode();
	int ndof = 2;              // up to now we only consider 2D problems!!!
	
	RCP<Epetra_SerialDenseMatrix> M_temp = rcp(new Epetra_SerialDenseMatrix(nrow,ncol));
	RCP<Epetra_SerialDenseMatrix> M_seg = rcp(new Epetra_SerialDenseMatrix(nrow*ndof,ncol*ndof));
	
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
		
		/* loop over all M_seg matrix entries
		   nrow represents the slave Lagrange multipliers !!!
		   ncol represents the master dofs !!!
		   (this DOES matter here for M_seg, as it might
		   sometimes be rectangular, not quadratic!)              */
		for (int j=0;j<nrow;++j)
		{
			for (int k=0;k<ncol;++k)
			{
				// multiply the two shape functions
				double prod = dualval[j]*mval[k];
				// add current Gauss point's contribution to M_seg  
				(*M_temp)(j,k) += prod*dxdsxi*dsxideta*wgt; 
			}
		}	
	} // for (int gp=0;gp<nGP();++gp)
	
	// fill M_seg matrix with M_temp matrix entries
	// (each M_temp value is multiplied with a (dof)-unit-matrix)
	for (int j=0;j<nrow*ndof;++j)
	{
		for (int k=0;k<ncol*ndof;++k)
		{
			int jindex = (int)(j/ndof);
			int kindex = (int)(k/ndof);
			// isolate the M_seg entries to be filled
			if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
				(*M_seg)(j,k) = (*M_temp)(jindex,kindex);
		}
	}

	return M_seg;
}

/*----------------------------------------------------------------------*
 |  Integrate a 1D slave / master overlap                     popp 01/08|
 |	This method integrates the modification to the Mortar matrix M      |
 |  for curved interface (Paper by Puso/Wohlmuth) from given local       |
 |	coordinates sxia to sxib. The corresponding master side local       |
 |  element coordinates given by mxia and mxib												  |
 |	Output is an Epetra_SerialDenseMatrix holding the int. values 			|
 *----------------------------------------------------------------------*/
RCP<Epetra_SerialDenseMatrix> CONTACT::Integrator::Integrate_Mmod(CONTACT::CElement& sele,
																					 										    double sxia, double sxib,
																					 										    CONTACT::CElement& mele,
																					 										    double mxia, double mxib)
{
	// check input data
	if ((!sele.IsSlave()) || (mele.IsSlave()))
		dserror("ERROR: Integrate_Mmod called on a wrong type of CElement pair!");
	if ((sxia<-1.0) || (sxib>1.0))
		dserror("ERROR: Integrate_Mmod called with infeasible slave limits!");
	if ((mxia<-1.0) || (mxib>1.0))
			dserror("ERROR: Integrate_Mmod called with infeasible master limits!");
	
	// create empty Mmod_seg object and wrap it with RCP
	int nrow  = sele.NumNode();
	int nrowdof = 2;							// up to now we only consider 2D problems!!!
	int ncol  = mele.NumNode();
	int ncoldof = 2;							// up to now we only consider 2D problems!!!
	
	RCP<Epetra_SerialDenseMatrix> Mmod_seg = rcp(new Epetra_SerialDenseMatrix(nrow*nrowdof,ncol*ncoldof));
	
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
			dserror("ERROR: Integrate_Mmod: Gauss point projection failed!");
		
		// evaluate trace space shape functions (on both elements)
		sele.EvaluateShape_1D(sxi,sval,sderiv,nrow);
		mele.EvaluateShape_1D(mxi,mval,mderiv,ncol);
		
		// build the delta function of slave side shape functions
		double delta_sval = sval[0]-sval[1];
		
		// evaluate the two slave side Jacobians
		double dxdsxi = sele.Jacobian_1D(sval,sderiv,scoord);
		double dsxideta = -0.5*sxia + 0.5*sxib;
		
		/* loop over all Mmod_seg matrix entries
		   nrow represents the slave Lagrange multipliers !!!
		   ncol represents the master dofs !!!
		   (this DOES matter here for Mmod_seg, as it might
		   sometimes be rectangular, not quadratic!)              */
		for (int j=0;j<nrow*nrowdof;++j)
		{
			for (int k=0;k<ncol*ncoldof;++k)
			{
				// multiply the two shape functions
				int mindex = (int)(k/ncoldof);
				double prod = 0.5*delta_sval*mval[mindex];
				// add current Gauss point's contribution to Mmod_seg  
				(*Mmod_seg)(j,k) += prod*dxdsxi*dsxideta*wgt; 
			}
		}	
	} // for (int gp=0;gp<nGP();++gp)
	
	/* OLD VERSION
	// prepare computation of purely geometric part of Mmod entries
	CNode* snode0 = static_cast<CNode*>(sele.Nodes()[0]);
	CNode* snode1 = static_cast<CNode*>(sele.Nodes()[1]);

	// normals
	double n[2][2];
	n[0][0] = snode0->n()[0];
	n[0][1] = snode0->n()[1];
	n[1][0] = snode1->n()[0];
	n[1][1] = snode1->n()[1];
	
	// tangents
	double t[2][2];
	t[0][0] = -n[0][1];
	t[0][1] =  n[0][0];
	t[1][0] = -n[1][1];
	t[1][1] =  n[1][0];
	
	// delta normals	
	double dn[2];
	dn[0] = n[0][0]-n[1][0];
	dn[1] = n[0][1]-n[1][1];
	
	// delta tangents	
	double dt[2];
	dt[0] = t[0][0]-t[1][0];
	dt[1] = t[0][1]-t[1][1];
	
	// loop over all Mmod_seg matrix entries
	for (int j=0;j<nrow*nrowdof;++j)
	{
		// prepare indices
		int snode = (int)(j/nrowdof);
		int sdof = (int)(j%nrowdof);
		
		for (int k=0;k<ncol*ncoldof;++k)
		{
			// prepare indices
			int mdof = (int)(k%ncoldof);
			
			// multiply geometric part onto Mmod  
			double val = n[snode][sdof] * dn[mdof] + t[snode][sdof] * dt[mdof];
			(*Mmod_seg)(j,k) *= val; 
		}
	}	
	*/ //OLD VERSION
	
	// NEW VERSION
	// prepare computation of purely geometric part of Mmod entries
	CNode* snode0 = static_cast<CNode*>(sele.Nodes()[0]);
	CNode* snode1 = static_cast<CNode*>(sele.Nodes()[1]);

	// normals
	double n[2][2];
	n[0][0] = snode0->n()[0];
	n[0][1] = snode0->n()[1];
	n[1][0] = snode1->n()[0];
	n[1][1] = snode1->n()[1];
	
	// tangents
	double t[2][2];
	t[0][0] = -n[0][1];
	t[0][1] =  n[0][0];
	t[1][0] = -n[1][1];
	t[1][1] =  n[1][0];
	
	// scalar product n1 * n2
	double n1n2 = 0.0;
	for (int i=0;i<3;++i)
		n1n2+=n[0][i]*n[1][i];
	
	// vector product n1 x n2
	double n1xn2 = n[0][0]*n[1][1] - n[0][1]*n[1][0];
	
	// // multiply geometric part onto Mmod  
	for (int i=0;i<ncol;++i)
	{
		(*Mmod_seg)(0,0+i*ncoldof) *=  (1.0-n1n2);
		(*Mmod_seg)(1,0+i*ncoldof) *=  n1xn2;
		(*Mmod_seg)(0,1+i*ncoldof) *= -n1xn2;
		(*Mmod_seg)(1,1+i*ncoldof) *=  (1.0-n1n2);
	
		(*Mmod_seg)(2,0+i*ncoldof) *=  (n1n2-1.0);
		(*Mmod_seg)(3,0+i*ncoldof) *=  n1xn2;
		(*Mmod_seg)(2,1+i*ncoldof) *= -n1xn2;
		(*Mmod_seg)(3,1+i*ncoldof) *=  (n1n2-1.0);
	}
	
	return Mmod_seg;
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
	
	// create empty g_seg object and wrap it with RCP
	int nrow = sele.NumNode();
	int ncol = mele.NumNode();
	RCP<Epetra_SerialDenseVector> g_seg = rcp(new Epetra_SerialDenseVector(nrow));
	
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
		if (length<1.0e-12) dserror("ERROR: Integrate_g: Divide by zero!");
		
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
		
		//if (gap<=0)
		//	gap = -sqrt(-gap);
		//else
		//	gap = sqrt(gap);
		
#ifdef DEBUG
		//cout << "GP gap: " << gap << endl;
#endif // #ifdef DEBUG
		
		// evaluate the two slave side Jacobians
		double dxdsxi = sele.Jacobian_1D(sval,sderiv,scoord);
		double dsxideta = -0.5*sxia + 0.5*sxib;
		
		/* loop over all g_seg vector entries
			 nrow represents the slave side dofs !!!                */
		for (int j=0;j<nrow;++j)
		{
			double prod = dualval[j]*gap;
			// add current Gauss point's contribution to g_seg  
			(*g_seg)(j) += prod*dxdsxi*dsxideta*wgt; 
		}
		
	} // for (int gp=0;gp<nGP();++gp)

	return g_seg;
}

/*----------------------------------------------------------------------*
 |  Assemble D contribution                                   popp 01/08|
 |	This method assembles the contrubution of a 1D slave element        |
 |  to the D map of the adjacent slave nodes.         			            |
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::Assemble_D(CONTACT::Interface& inter,
																		 CONTACT::CElement& sele,
																		 Epetra_SerialDenseMatrix& D_seg)
{
	/*
#ifdef DEBUG
	cout << "Calling proc: " << inter.Comm().MyPID() << endl;
	cout << D_seg << endl;
#endif // #ifdef DEBUG
	*/
	// get adjacent nodes to assemble to
	DRT::Node** snodes = sele.Nodes();
	if (!snodes)
		dserror("ERROR: Assemble_D: Null pointer for snodes!");
	
	// loop over all slave nodes
	for (int slave=0;slave<sele.NumNode();++slave)
	{
		CONTACT::CNode* snode = static_cast<CONTACT::CNode*>(snodes[slave]);
		//const int* sdofs = snode->Dofs();
		int sndof = snode->NumDof();
		
		// only process slave node rows that belong to this proc
		if (snode->Owner() != inter.Comm().MyPID())
			continue;
		
		// loop over all dofs of the slave node
		for (int sdof=0;sdof<sndof;++sdof)
		{
			// loop over all slave nodes again ("master nodes")
			for (int master=0;master<sele.NumNode();++master)
			{
				CONTACT::CNode* mnode = static_cast<CONTACT::CNode*>(snodes[master]);
				const int* mdofs = mnode->Dofs();
				int mndof = mnode->NumDof();
				
				// loop over all dofs of the slave node again ("master dofs")
				for (int mdof=0;mdof<mndof;++mdof)
				{
					int col = mdofs[mdof];
					double val = D_seg(slave*sndof+sdof,master*mndof+mdof);
					snode->AddDValue(sdof,col,val);
				}
			}
		}
		/*
#ifdef DEBUG
		cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << endl;
		map<int, double> nodemap0 = (snode->GetD())[0];
		map<int, double> nodemap1 = (snode->GetD())[1];
		typedef map<int,double>::const_iterator CI;
				
		cout << "Row dof id: " << sdofs[0] << endl;;
		for (CI p=nodemap0.begin();p!=nodemap0.end();++p)
			cout << p->first << '\t' << p->second << endl;
				
		cout << "Row dof id: " << sdofs[1] << endl;
		for (CI p=nodemap1.begin();p!=nodemap1.end();++p)
			cout << p->first << '\t' << p->second << endl;
#endif // #ifdef DEBUG
		*/
	}
	
	return true;
}

/*----------------------------------------------------------------------*
 |  Assemble M contribution                                   popp 01/08|
 |	This method assembles the contrubution of a 1D slave / master			  |
 |  overlap pair to the M map of the adjacent slave nodes.         			|
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::Assemble_M(CONTACT::Interface& inter,
																		 CONTACT::CElement& sele,
																		 CONTACT::CElement& mele,
																		 Epetra_SerialDenseMatrix& M_seg)
{
	/*
#ifdef DEBUG
	cout << "Calling proc: " << inter.Comm().MyPID() << endl;
	cout << M_seg << endl;
#endif // #ifdef DEBUG
	*/
	// get adjacent slave nodes and master nodes
	DRT::Node** snodes = sele.Nodes();
	if (!snodes)
		dserror("ERROR: Assemble_M: Null pointer for snodes!");
	DRT::Node** mnodes = mele.Nodes();
	if (!mnodes)
		dserror("ERROR: Assemble_M: Null pointer for mnodes!");
	
	// loop over all slave nodes
	for (int slave=0;slave<sele.NumNode();++slave)
	{
		CONTACT::CNode* snode = static_cast<CONTACT::CNode*>(snodes[slave]);
		//const int* sdofs = snode->Dofs();
		int sndof = snode->NumDof();
		
		// only process slave node rows that belong to this proc
		if (snode->Owner() != inter.Comm().MyPID())
			continue;
		
		// loop over all dofs of the slave node
		for (int sdof=0;sdof<sndof;++sdof)
		{
			// loop over all master nodes
			for (int master=0;master<mele.NumNode();++master)
			{
				CONTACT::CNode* mnode = static_cast<CONTACT::CNode*>(mnodes[master]);
				const int* mdofs = mnode->Dofs();
				int mndof = mnode->NumDof();
				
				// loop over all dofs of the master node
				for (int mdof=0;mdof<mndof;++mdof)
				{
					int col = mdofs[mdof];
					double val = M_seg(slave*sndof+sdof,master*mndof+mdof);
					snode->AddMValue(sdof,col,val);
				}
			}
		}
		/*
#ifdef DEBUG
		cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << endl;
		map<int, double> nodemap0 = (snode->GetM())[0];
		map<int, double> nodemap1 = (snode->GetM())[1];
		typedef map<int,double>::const_iterator CI;
		
		cout << "Row dof id: " << sdofs[0] << endl;;
		for (CI p=nodemap0.begin();p!=nodemap0.end();++p)
			cout << p->first << '\t' << p->second << endl;
		
		cout << "Row dof id: " << sdofs[1] << endl;
		for (CI p=nodemap1.begin();p!=nodemap1.end();++p)
			cout << p->first << '\t' << p->second << endl;
#endif // #ifdef DEBUG
		 */
	}
	
	return true;
}

/*----------------------------------------------------------------------*
 |  Assemble Mmod contribution                                popp 01/08|
 |	This method assembles the contribution of a 1D slave / master			  |
 |  overlap pair to the Mmod map of the adjacent slave nodes.      			|
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::Assemble_Mmod(CONTACT::Interface& inter,
																		    CONTACT::CElement& sele,
																		    CONTACT::CElement& mele,
																		    Epetra_SerialDenseMatrix& Mmod_seg)
{
	/*
#ifdef DEBUG
	cout << "Calling proc: " << inter.Comm().MyPID() << endl;
	cout << Mmod_seg << endl;
#endif // #ifdef DEBUG
	*/
	
	// get adjacent slave nodes and master nodes
	DRT::Node** snodes = sele.Nodes();
	if (!snodes)
		dserror("ERROR: Assemble_Mmod: Null pointer for snodes!");
	DRT::Node** mnodes = mele.Nodes();
	if (!mnodes)
		dserror("ERROR: Assemble_Mmod: Null pointer for mnodes!");
	
	// loop over all slave nodes
	for (int slave=0;slave<sele.NumNode();++slave)
	{
		CONTACT::CNode* snode = static_cast<CONTACT::CNode*>(snodes[slave]);
		//const int* sdofs = snode->Dofs();
		int sndof = snode->NumDof();
		
		// only process slave node rows that belong to this proc
		if (snode->Owner() != inter.Comm().MyPID())
			continue;
		
		// loop over all dofs of the slave node
		for (int sdof=0;sdof<sndof;++sdof)
		{
			// loop over all master nodes
			for (int master=0;master<mele.NumNode();++master)
			{
				CONTACT::CNode* mnode = static_cast<CONTACT::CNode*>(mnodes[master]);
				const int* mdofs = mnode->Dofs();
				int mndof = mnode->NumDof();
				
				// loop over all dofs of the master node
				for (int mdof=0;mdof<mndof;++mdof)
				{
					int col = mdofs[mdof];
					double val = Mmod_seg(slave*sndof+sdof,master*mndof+mdof);
					snode->AddMmodValue(sdof,col,val);
				}
			}
		}
		/*
#ifdef DEBUG
		cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << endl;
		map<int, double> nodemap0 = (snode->GetMmod())[0];
		map<int, double> nodemap1 = (snode->GetMmod())[1];
		typedef map<int,double>::const_iterator CI;
		
		cout << "Row dof id: " << sdofs[0] << endl;;
		for (CI p=nodemap0.begin();p!=nodemap0.end();++p)
			cout << p->first << '\t' << p->second << endl;
		
		cout << "Row dof id: " << sdofs[1] << endl;
		for (CI p=nodemap1.begin();p!=nodemap1.end();++p)
			cout << p->first << '\t' << p->second << endl;
#endif // #ifdef DEBUG
		 */
	}
	
	return true;
}

/*----------------------------------------------------------------------*
 |  Assemble g~ contribution                                  popp 01/08|
 |	This method assembles the contribution of a 1D slave / master			  |
 |  overlap pair to the weighted gap of the adjacent slave nodes.     	|
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::Assemble_g(CONTACT::Interface& inter,
																		 CONTACT::CElement& sele,
																		 Epetra_SerialDenseVector& g_seg)
{
	/*
#ifdef DEBUG
	cout << "Calling proc: " << inter.Comm().MyPID() << endl;
	cout << g_seg << endl;
#endif // #ifdef DEBUG
	*/
	
	// get adjacent slave to assemble to
	DRT::Node** snodes = sele.Nodes();
	if (!snodes)
		dserror("ERROR: Assemble_g: Null pointer for snodes!");
	
	// loop over all slave nodes
	for (int slave=0;slave<sele.NumNode();++slave)
	{
		CONTACT::CNode* snode = static_cast<CONTACT::CNode*>(snodes[slave]);
		
		// only process slave node rows that belong to this proc
		if (snode->Owner() != inter.Comm().MyPID())
			continue;
		
		double val = g_seg(slave);
		snode->AddgValue(val);
		/*
#ifdef DEBUG
		cout << "Node: " << snode->Id() << "  Owner: " << snode->Owner() << endl;
		cout << "Weighted gap: " << snode->Getg() << endl;
#endif // #ifdef DEBUG
		*/
	}
	
	return true;
}

#endif //#ifdef CCADISCRET
