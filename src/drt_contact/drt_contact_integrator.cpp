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
CONTACT::Integrator::Integrator(int ngp, bool oned) :
oned_(oned),
ngp_(ngp)
{
	if (oned)
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
		
		
	} // if (oned)
	
	else
		dserror("ERROR: Integrator: 2D case not yet implemented!");
}

/*----------------------------------------------------------------------*
 |  Integrate a 1D slave element overlap                      popp 01/08|
 |	This method integrates 2 functions on the same (slave) CEelement    |
 |	from given local coordinates sxia to sxib														|
 |	Output is an Epetra_SerialDenseMatrix holding the int. values 			|
 *----------------------------------------------------------------------*/
RCP<Epetra_SerialDenseMatrix> CONTACT::Integrator::IntegrateD(CONTACT::CElement& sele,
																					 										double sxia, double sxib)
{
	//check input data
	if (!sele.IsSlave())
		dserror("ERROR: IntegrateD called on a non-slave CElement!");
	if ((sxia<-1.0) || (sxib>1.0))
		dserror("ERROR: IntegrateD called with infeasible slave limits!");
	
	// create empty dseg object and wrap it with RCP
	int nrow = sele.NumNode();
	int ndof = 2;              // up to now we only consider 2D problems!!!
	int ncol = nrow;
	
	RCP<Epetra_SerialDenseMatrix> dtemp = rcp(new Epetra_SerialDenseMatrix(nrow,ncol));
	RCP<Epetra_SerialDenseMatrix> dseg = rcp(new Epetra_SerialDenseMatrix(nrow*ndof,ncol*ndof));
	
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
		sele.EvaluateShape1D(sxi,val,deriv,nrow);
		sele.EvaluateShapeDual1D(sxi,dualval,dualderiv,nrow);
		
		// evaluate the two Jacobians
		double dxdsxi = sele.Jacobian1D(val,deriv,coord);
		double dsxideta = -0.5*sxia + 0.5*sxib;
		
		/* loop over all dtemp matrix entries
		   nrow represents the dofs !!!
		   ncol represents the Lagrange multipliers !!!
		   (although this does not really matter here for dseg,
		   as it will turn out to be diagonal anyway)             */
		for (int j=0;j<nrow;++j)
		{
			for (int k=0;k<ncol;++k)
			{
				// multiply the two shape functions
				double prod = val[j]*dualval[k];
				// add current Gauss point's contribution to dtemp  
				(*dtemp)(j,k) += prod*dxdsxi*dsxideta*wgt; 
			}
		}	
	} // for (int gp=0;gp<nGP();++gp)

	// fill dseg matrix with dtemp matrix entries
	// (each dtemp value is multiplied with a (dof)-unit-matrix)
	for (int j=0;j<nrow*ndof;++j)
	{
		for (int k=0;k<ncol*ndof;++k)
		{
			int jindex = (int)(j/ndof);
			int kindex = (int)(k/ndof);
			// isolate the dseg entries to be filled
			if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
				(*dseg)(j,k) = (*dtemp)(jindex,kindex);
		}
	}

	return dseg;
}

/*----------------------------------------------------------------------*
 |  Integrate a 1D slave / master overlap                     popp 01/08|
 |	This method integrates a slave side function (dual shape fct.)      |
 |  and a master side function (standard shape fct.) from given local   |
 |	coordinates sxia to sxib. The corresponding master side local       |
 |  element coordinates given by mxia and mxib												  |
 |	Output is an Epetra_SerialDenseMatrix holding the int. values 			|
 *----------------------------------------------------------------------*/
RCP<Epetra_SerialDenseMatrix> CONTACT::Integrator::IntegrateM(CONTACT::CElement& sele,
																					 										double sxia, double sxib,
																					 										CONTACT::CElement& mele,
																					 										double mxia, double mxib)
{
	// check input data
	if ((!sele.IsSlave()) || (mele.IsSlave()))
		dserror("ERROR: IntegrateM called on a wrong type of CElement pair!");
	if ((sxia<-1.0) || (sxib>1.0))
		dserror("ERROR: IntegrateM called with infeasible slave limits!");
	if ((mxia<-1.0) || (mxib>1.0))
			dserror("ERROR: IntegrateM called with infeasible master limits!");
	
	// create empty mseg object and wrap it with RCP
	int nrow = sele.NumNode();
	int ncol = mele.NumNode();
	int ndof = 2;              // up to now we only consider 2D problems!!!
	
	RCP<Epetra_SerialDenseMatrix> mtemp = rcp(new Epetra_SerialDenseMatrix(nrow,ncol));
	RCP<Epetra_SerialDenseMatrix> mseg = rcp(new Epetra_SerialDenseMatrix(nrow*ndof,ncol*ndof));
	
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
		projector.ProjectGaussPoint(sele,sxi,mele,mxi);
		
		// check GP projection
		if ((mxi[0]<mxia) || (mxi[0]>mxib))
		{
			cout << "Slave ID: " << sele.Id() << " Master ID: " << mele.Id() << endl;
			cout << "Slave nodes: " << sele.NodeIds()[0] << " " << sele.NodeIds()[1] << endl;
			cout << "Master nodes: " << mele.NodeIds()[0] << " " << mele.NodeIds()[1] << endl;
			cout << "sxia: " << sxia << " sxib: " << sxib << endl;
			cout << "mxia: " << mxia << " mxib: " << mxib << endl;
			dserror("ERROR: IntegrateM: Gauss point projection failed! mxi=%d",mxi[0]);
		}
		// evaluate dual space shape functions (on slave element)
		sele.EvaluateShapeDual1D(sxi,dualval,dualderiv,nrow);
		
		// evaluate trace space shape functions (on both elements)
		sele.EvaluateShape1D(sxi,sval,sderiv,nrow);
		mele.EvaluateShape1D(mxi,mval,mderiv,ncol);
		
		// evaluate the two slave side Jacobians
		double dxdsxi = sele.Jacobian1D(sval,sderiv,scoord);
		double dsxideta = -0.5*sxia + 0.5*sxib;
		
		/* loop over all mseg matrix entries
		   nrow represents the slave Lagrange multipliers !!!
		   ncol represents the master dofs !!!
		   (this DOES matter here for mseg, as it might
		   sometimes be rectangular, not quadratic!)              */
		for (int j=0;j<nrow;++j)
		{
			for (int k=0;k<ncol;++k)
			{
				// multiply the two shape functions
				double prod = dualval[j]*mval[k];
				// add current Gauss point's contribution to mseg  
				(*mtemp)(j,k) += prod*dxdsxi*dsxideta*wgt; 
			}
		}	
	} // for (int gp=0;gp<nGP();++gp)
	
	// fill mseg matrix with mtemp matrix entries
	// (each mtemp value is multiplied with a (dof)-unit-matrix)
	for (int j=0;j<nrow*ndof;++j)
	{
		for (int k=0;k<ncol*ndof;++k)
		{
			int jindex = (int)(j/ndof);
			int kindex = (int)(k/ndof);
			// isolate the mseg entries to be filled
			if ((j==k) || ((j-jindex*ndof)==(k-kindex*ndof)))
				(*mseg)(j,k) = (*mtemp)(jindex,kindex);
		}
	}

	return mseg;
}

/*----------------------------------------------------------------------*
 |  Integrate a 1D slave / master overlap                     popp 01/08|
 |	This method integrates the modification to the Mortar matrix M      |
 |  for curved interface (Paper by Puso/Wohlmuth) from given local       |
 |	coordinates sxia to sxib. The corresponding master side local       |
 |  element coordinates given by mxia and mxib												  |
 |	Output is an Epetra_SerialDenseMatrix holding the int. values 			|
 *----------------------------------------------------------------------*/
RCP<Epetra_SerialDenseMatrix> CONTACT::Integrator::IntegrateMmod(CONTACT::CElement& sele,
																					 										    double sxia, double sxib,
																					 										    CONTACT::CElement& mele,
																					 										    double mxia, double mxib)
{
	// check input data
	if ((!sele.IsSlave()) || (mele.IsSlave()))
		dserror("ERROR: IntegrateMmod called on a wrong type of CElement pair!");
	if ((sxia<-1.0) || (sxib>1.0))
		dserror("ERROR: IntegrateMmod called with infeasible slave limits!");
	if ((mxia<-1.0) || (mxib>1.0))
			dserror("ERROR: IntegrateMmod called with infeasible master limits!");
	
	// create empty mmodseg object and wrap it with RCP
	int nrow  = sele.NumNode();
	int nrowdof = 2;							// up to now we only consider 2D problems!!!
	int ncol  = mele.NumNode();
	int ncoldof = 2;							// up to now we only consider 2D problems!!!
	
	RCP<Epetra_SerialDenseMatrix> mmodseg = rcp(new Epetra_SerialDenseMatrix(nrow*nrowdof,ncol*ncoldof));
	
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
		projector.ProjectGaussPoint(sele,sxi,mele,mxi);
		
		// check GP projection
		if ((mxi[0]<mxia) || (mxi[0]>mxib))
			dserror("ERROR: IntegrateMmod: Gauss point projection failed!");
		
		// evaluate trace space shape functions (on both elements)
		sele.EvaluateShape1D(sxi,sval,sderiv,nrow);
		mele.EvaluateShape1D(mxi,mval,mderiv,ncol);
		
		// build the delta function of slave side shape functions
		double deltasval = sval[0]-sval[1];
		
		// evaluate the two slave side Jacobians
		double dxdsxi = sele.Jacobian1D(sval,sderiv,scoord);
		double dsxideta = -0.5*sxia + 0.5*sxib;
		
		/* loop over all mmodseg matrix entries
		   nrow represents the slave Lagrange multipliers !!!
		   ncol represents the master dofs !!!
		   (this DOES matter here for mmodseg, as it might
		   sometimes be rectangular, not quadratic!)              */
		for (int j=0;j<nrow*nrowdof;++j)
		{
			for (int k=0;k<ncol*ncoldof;++k)
			{
				// multiply the two shape functions
				int mindex = (int)(k/ncoldof);
				double prod = 0.5*deltasval*mval[mindex];
				// add current Gauss point's contribution to mmodseg  
				(*mmodseg)(j,k) += prod*dxdsxi*dsxideta*wgt; 
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
	
	// loop over all mmodseg matrix entries
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
			(*mmodseg)(j,k) *= val; 
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
		(*mmodseg)(0,0+i*ncoldof) *=  (1.0-n1n2);
		(*mmodseg)(1,0+i*ncoldof) *=  n1xn2;
		(*mmodseg)(0,1+i*ncoldof) *= -n1xn2;
		(*mmodseg)(1,1+i*ncoldof) *=  (1.0-n1n2);
	
		(*mmodseg)(2,0+i*ncoldof) *=  (n1n2-1.0);
		(*mmodseg)(3,0+i*ncoldof) *=  n1xn2;
		(*mmodseg)(2,1+i*ncoldof) *= -n1xn2;
		(*mmodseg)(3,1+i*ncoldof) *=  (n1n2-1.0);
	}
	
	return mmodseg;
}

/*----------------------------------------------------------------------*
 |  Integrate gap on a 1D slave / master overlap              popp 01/08|
 |	This method integrates a slave side function (dual shape fct.)      |
 |  and the gap function g = ( ( sx - mx ) * n )  from given local      |
 |	coordinates sxia to sxib. The corresponding master side local       |
 |  element coordinates given by mxia and mxib												  |
 |	Output is an Epetra_SerialDenseVector holding the int. values 			|
 *----------------------------------------------------------------------*/
RCP<Epetra_SerialDenseVector> CONTACT::Integrator::IntegrateG(CONTACT::CElement& sele,
																					 										double sxia, double sxib,
																					 										CONTACT::CElement& mele,
																					 										double mxia, double mxib)
{
	// check input data
	if ((!sele.IsSlave()) || (mele.IsSlave()))
		dserror("ERROR: IntegrateG called on a wrong type of CElement pair!");
	if ((sxia<-1.0) || (sxib>1.0))
		dserror("ERROR: IntegrateG called with infeasible slave limits!");
	if ((mxia<-1.0) || (mxib>1.0))
			dserror("ERROR: IntegrateG called with infeasible master limits!");
	
	// create empty gseg object and wrap it with RCP
	int nrow = sele.NumNode();
	int ncol = mele.NumNode();
	RCP<Epetra_SerialDenseVector> gseg = rcp(new Epetra_SerialDenseVector(nrow));
	
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
	if(!mynodes) dserror("ERROR: IntegrateG: Null pointer!");
	
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
		projector.ProjectGaussPoint(sele,sxi,mele,mxi);
		
		// check GP projection
		if ((mxi[0]<mxia) || (mxi[0]>mxib))
			dserror("ERROR: IntegrateG: Gauss point projection failed!");
		
		// evaluate dual space shape functions (on slave element)
		sele.EvaluateShapeDual1D(sxi,dualval,dualderiv,nrow);
		
		// evaluate trace space shape functions (on both elements)
		sele.EvaluateShape1D(sxi,sval,sderiv,nrow);
		mele.EvaluateShape1D(mxi,mval,mderiv,ncol);
		
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
		if (length<1.0e-12) dserror("ERROR: IntegrateG: Divide by zero!");
		
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
		double dxdsxi = sele.Jacobian1D(sval,sderiv,scoord);
		double dsxideta = -0.5*sxia + 0.5*sxib;
		
		/* loop over all gseg vector entries
			 nrow represents the slave side dofs !!!                */
		for (int j=0;j<nrow;++j)
		{
			double prod = dualval[j]*gap;
			// add current Gauss point's contribution to gseg  
			(*gseg)(j) += prod*dxdsxi*dsxideta*wgt; 
		}
		
	} // for (int gp=0;gp<nGP();++gp)

	return gseg;
}

/*----------------------------------------------------------------------*
 |  Assemble D contribution                                   popp 01/08|
 |	This method assembles the contrubution of a 1D slave element        |
 |  to the D map of the adjacent slave nodes.         			            |
 *----------------------------------------------------------------------*/
bool CONTACT::Integrator::AssembleD(CONTACT::Interface& inter,
																		 CONTACT::CElement& sele,
																		 Epetra_SerialDenseMatrix& dseg)
{
	/*
#ifdef DEBUG
	cout << "Calling proc: " << inter.Comm().MyPID() << endl;
	cout << dseg << endl;
#endif // #ifdef DEBUG
	*/
	// get adjacent nodes to assemble to
	DRT::Node** snodes = sele.Nodes();
	if (!snodes)
		dserror("ERROR: AssembleD: Null pointer for snodes!");
	
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
					double val = dseg(slave*sndof+sdof,master*mndof+mdof);
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
bool CONTACT::Integrator::AssembleM(CONTACT::Interface& inter,
																		 CONTACT::CElement& sele,
																		 CONTACT::CElement& mele,
																		 Epetra_SerialDenseMatrix& mseg)
{
	/*
#ifdef DEBUG
	cout << "Calling proc: " << inter.Comm().MyPID() << endl;
	cout << mseg << endl;
#endif // #ifdef DEBUG
	*/
	// get adjacent slave nodes and master nodes
	DRT::Node** snodes = sele.Nodes();
	if (!snodes)
		dserror("ERROR: AssembleM: Null pointer for snodes!");
	DRT::Node** mnodes = mele.Nodes();
	if (!mnodes)
		dserror("ERROR: AssembleM: Null pointer for mnodes!");
	
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
					double val = mseg(slave*sndof+sdof,master*mndof+mdof);
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
bool CONTACT::Integrator::AssembleMmod(CONTACT::Interface& inter,
																		    CONTACT::CElement& sele,
																		    CONTACT::CElement& mele,
																		    Epetra_SerialDenseMatrix& mmodseg)
{
	/*
#ifdef DEBUG
	cout << "Calling proc: " << inter.Comm().MyPID() << endl;
	cout << mmodseg << endl;
#endif // #ifdef DEBUG
	*/
	
	// get adjacent slave nodes and master nodes
	DRT::Node** snodes = sele.Nodes();
	if (!snodes)
		dserror("ERROR: AssembleMmod: Null pointer for snodes!");
	DRT::Node** mnodes = mele.Nodes();
	if (!mnodes)
		dserror("ERROR: AssembleMmod: Null pointer for mnodes!");
	
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
					double val = mmodseg(slave*sndof+sdof,master*mndof+mdof);
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
bool CONTACT::Integrator::AssembleG(CONTACT::Interface& inter,
																		 CONTACT::CElement& sele,
																		 Epetra_SerialDenseVector& gseg)
{
	/*
#ifdef DEBUG
	cout << "Calling proc: " << inter.Comm().MyPID() << endl;
	cout << gseg << endl;
#endif // #ifdef DEBUG
	*/
	
	// get adjacent slave to assemble to
	DRT::Node** snodes = sele.Nodes();
	if (!snodes)
		dserror("ERROR: AssembleG: Null pointer for snodes!");
	
	// loop over all slave nodes
	for (int slave=0;slave<sele.NumNode();++slave)
	{
		CONTACT::CNode* snode = static_cast<CONTACT::CNode*>(snodes[slave]);
		
		// only process slave node rows that belong to this proc
		if (snode->Owner() != inter.Comm().MyPID())
			continue;
		
		double val = gseg(slave);
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
