/*!----------------------------------------------------------------------
\file fluid3_xfem_evaluate.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "fluid3_xfem.H"
#include "fluid3_xfem_enrichment.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"


using namespace Enrichments;
using namespace DRT::Utils;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::XFluid3::Evaluate(ParameterList&           params,
                                     DRT::Discretization&      discretization,
                                     vector<int>&              lm,
                                     Epetra_SerialDenseMatrix& elemat1,
                                     Epetra_SerialDenseMatrix& elemat2,
                                     Epetra_SerialDenseVector& elevec1,
                                     Epetra_SerialDenseVector& elevec2,
                                     Epetra_SerialDenseVector& elevec3)
{
    // get the action required
    const string action = params.get<string>("action","none");
    
    const DRT::Elements::XFluid3::ActionType act = convertStringToActionType(action);

    // get the material
    const MATERIAL* actmat = &(mat[material_-1]);

  	switch(act)
  	{
   	case calc_fluid_systemmat_and_residual:
    {
        const int numenr = 1;
        const int numnode = NumNode();
        
        // define once and for all the index mapping within this element call
        vector<ElementEnrichment>  enrichments;
        
        //loop over nodes to find all enrichment id's
        // for now only standard enrichment shall be available
        ElementEnrichment ele_enr = ElementEnrichment(0, enrichment_type_standard);
        enrichments.push_back(ele_enr);
        
        map<drtIdx, int> dofmap;
        
        vector<DRT::Node*> nodes;
        for (int inode=0;inode<numnode;++inode)
        {
            DRT::Node& noderef = *(Nodes()[inode]);
            DRT::Node* node_ptr = Nodes()[inode];
//            noderef.Print(cout);
            
            //cout << &(node.Print());
            //cout << (&node_ptr->Print());
        }
//        this->Print(cout);
//        vector<DRT::Node*>::const_iterator curr;
//        for (curr=nodes.begin(); curr!=nodes.end(); ++curr)
//            cout << "      " << &(curr) << endl;
        
        //vector<DRT::Node> nodes = Nodes();
//        for (int inode=0;inode<numnode;++inode)
//        {
//           const DRT::Node node = 
//        
//        for (int ienr=0;ienr<numenr;++ienr)
//        {
//            const ElementEnrichment eleenr = enrichments[ienr];
//            //cout << eleenr.getId() << endl;
//
////                if
//
//
//        }
//        
//        //                for (int idf=0;idf<NDF_;++idf)
////                {
////                    const drtIdx idx(ienr, inode, idf);
        

        //dserror("yo");


        
        // need current velocity and history vector
        RefCountPtr<const Epetra_Vector> vel_pre_np = discretization.GetState("u and p at time n+1 (trial)");
        dsassert(vel_pre_np!=null, "Cannot get state vectors 'velnp'");
        RefCountPtr<const Epetra_Vector> hist  = discretization.GetState("old solution data for rhs");
        dsassert(hist!=null, "Cannot get state vectors 'hist'");
      
        // extract local values from the global vectors
        vector<double> my_vel_pre_np(lm.size());
        DRT::Utils::ExtractMyValues(*vel_pre_np,my_vel_pre_np,lm);
        vector<double> myhist(lm.size());
        DRT::Utils::ExtractMyValues(*hist,myhist,lm);

        // split "my_vel_pre_np" into velocity part "myvelnp" and pressure part "myprenp"
        // Additionally only the velocity components of myhist are important!
        
        vector<double> myprenp(numnode);
        vector<double> myvelnp(NSD_*numnode);
        vector<double> myvhist(NSD_*numnode);
        //EnrEleVector myvhist(numenr);
      
        for (int i=0;i<numenr;++i)
        {
            
            for (int i=0;i<numnode;++i)
            {
                myvelnp[0+(i*NSD_)]=my_vel_pre_np[0+(i*4)];
                myvelnp[1+(i*NSD_)]=my_vel_pre_np[1+(i*4)];
                myvelnp[2+(i*NSD_)]=my_vel_pre_np[2+(i*4)];

                myprenp[i]=my_vel_pre_np[NSD_+(i*4)];

                myvhist[0+(i*NSD_)]=myhist[0+(i*4)];
                myvhist[1+(i*NSD_)]=myhist[1+(i*4)];
                myvhist[2+(i*NSD_)]=myhist[2+(i*4)];
                //myvhist[0][0][i] = myhist[0+(i*4)];
                //myvhist[0][1][i] = myhist[1+(i*4)];
                //myvhist[0][2][i] = myhist[2+(i*4)];
            }
        }

        // calculate element coefficient matrix and rhs       
        f3_sys_mat(lm,myvelnp,myprenp,myvhist,&elemat1,&elevec1,actmat,params);
        break;
    }
    case calc_fluid_beltrami_error:
    {
        // add error only for elements which are not ghosted
        if(this->Owner() == discretization.Comm().MyPID())
        {
        
          	// need current velocity and history vector
          	RefCountPtr<const Epetra_Vector> vel_pre_np = discretization.GetState("u and p at time n+1 (converged)");
            dsassert(vel_pre_np!=null, "Cannot get state vectors 'velnp'");
      
          	// extract local values from the global vectors
          	vector<double> my_vel_pre_np(lm.size());
          	DRT::Utils::ExtractMyValues(*vel_pre_np,my_vel_pre_np,lm);

          	// split "my_vel_pre_np" into velocity part "myvelnp" and pressure part "myprenp"
          	int numnode = NumNode();
          	vector<double> myprenp(numnode);
          	vector<double> myvelnp(NSD_*numnode);
      
          	for (int i=0;i<numnode;++i)
          	{
            	myvelnp[0+(i*NSD_)]=my_vel_pre_np[0+(i*4)];
            	myvelnp[1+(i*NSD_)]=my_vel_pre_np[1+(i*4)];
            	myvelnp[2+(i*NSD_)]=my_vel_pre_np[2+(i*4)];

            	myprenp[i]=my_vel_pre_np[NSD_+(i*4)];
          	}

          	// integrate beltrami error
          	f3_int_beltrami_err(myvelnp,myprenp,actmat,params);
        }
        break;
    }
    case calc_Shapefunction:
    {
     	shape_function_3D(elevec1,elevec2[0],elevec2[1],elevec2[2],this->Shape());
      	break;
    }
    case calc_ShapeDeriv1:
    {
    	shape_function_3D_deriv1(elemat1,elevec2[0],elevec2[1],elevec2[2],this->Shape());
     	break;
    }
    case calc_ShapeDeriv2:
    {
		shape_function_3D_deriv2(elemat2,elevec2[0],elevec2[1],elevec2[2],this->Shape());
      	break;
    }
    default:
        dserror("Unknown type of action for Fluid3");
    } // end of switch(act)

    return 0;
} // end of DRT::Elements::Fluid3::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      gammi 04/07|
 |                                                                      |
 |  The function is just a dummy. For the fluid elements, the           |
 |  integration of the volume neumann loads takes place in the element. |
 |  We need it there for the stabilisation terms!                       |
 *----------------------------------------------------------------------*/
int DRT::Elements::XFluid3::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


/*----------------------------------------------------------------------*
  |  calculate system matrix and rhs (private)                g.bau 03/07|
  *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3::f3_sys_mat(const vector<int>&        lm,
                                        const vector<double>&     evelnp,
                                        const vector<double>&     eprenp,
                                        const vector<double>&     evhist,
                                        Epetra_SerialDenseMatrix* sys_mat,
                                        Epetra_SerialDenseVector* residual,
                                        const struct _MATERIAL*   material,
                                        ParameterList&            params)
{
    const int numnode = this->NumNode();
    const DiscretizationType distype = this->Shape();
    
    const int numenr = 1; // number of enrichments within element (will be avariable later on)
    const int NDF = 4; // number of physical dof
    double  mat_ele[numenr][NDF][numnode][numenr][NDF][numnode]; // system matrix that will fill the epetra matrix
    double  vec_ele[numenr][NDF][numnode]; // system vector (residual) that will fill the epetra matrix

    
    // Fluid3 element using USFEM
    // get element node coordinates
    Epetra_SerialDenseMatrix xyze(NSD_, numnode);
    for(int inode=0;inode<numnode;inode++)
    {
        xyze(0,inode)=Nodes()[inode]->X()[0];
        xyze(1,inode)=Nodes()[inode]->X()[1];
        xyze(2,inode)=Nodes()[inode]->X()[2];
    }

    // dead load in element nodes
    Epetra_SerialDenseMatrix bodyforce(NSD_,numnode);
    this->f3_getbodyforce(bodyforce,numnode,params);

    /*---------------------------------------------- get viscosity ---*/
    // check here, if we really have a fluid !!
    dsassert(material->mattyp == m_fluid, "Material law is not of type m_fluid.");
    const double  visc = material->m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    // USFEM stabilization is default. No switch here at the moment.

    /*----------------------------------------- declaration of variables ---*/
    Epetra_SerialDenseVector    funct(numnode);
    Epetra_SerialDenseMatrix    deriv(NSD_,numnode);
    Epetra_SerialDenseMatrix    deriv2(6,numnode);
    double  enrshp[12][numenr][NDF][numnode];   //enriched shape function
        
    Epetra_SerialDenseMatrix    vderxy(NSD_,NSD_);
    vector<double>              pderxy(NSD_);
    Epetra_SerialDenseMatrix    vderxy2(NSD_,6);
    Epetra_SerialDenseMatrix    derxy(NSD_,numnode);
    Epetra_SerialDenseMatrix    derxy2(6,numnode);
    // data at integration point
    vector<double>              edeadng(NSD_); //
    vector<double>              histvec(NSD_); // history values
    vector<double>              gradp(NSD_);   // pressure gradient
    vector<double>              velint(NSD_);  // velocity
    vector<double>              gridvelint(NSD_); // grid velocity

    // get control parameter
    const double timefac=params.get<double>("time constant for integration",0.0);

    const bool is_stationary = params.get<bool>("using stationary formulation",false);
    
    // stabilization parameter
    const vector<double> tau = f3_caltau(xyze, evelnp, distype, visc, numnode, timefac, is_stationary);

    // flag for higher order elements
    const bool higher_order_ele = is_higher_order_element(distype);
    
    // gaussian points
    const IntegrationPoints3D  intpoints = getIntegrationPoints3D(gaussrule_);
        
    // calculate element volume to check gauss integration
    double  element_volume = 0.0;
    // start loop over integration points
    for (int iquad=0;iquad<intpoints.nquad;iquad++)
    {
        const double e1 = intpoints.qxg[iquad][0];
        const double e2 = intpoints.qxg[iquad][1];
        const double e3 = intpoints.qxg[iquad][2];
        shape_function_3D(funct,e1,e2,e3,distype);
        shape_function_3D_deriv1(deriv,e1,e2,e3,distype);
        if (higher_order_ele)
        {  
            shape_function_3D_deriv2(deriv2,e1,e2,e3,distype);
        }
    
        // compute Jacobian matrix and its determinante
        static const Epetra_SerialDenseMatrix xjm = getJacobiMatrix(xyze,deriv,numnode);
        const double det = getDeterminante(xjm);

        // compute first global derivates of shape functions
        f3_gder(derxy,deriv,xjm,det,numnode);

        // compute second global derivative of shape functions
        if (higher_order_ele)
        {
            f3_gder2(xyze,xjm,derxy,derxy2,deriv2,numnode);
        }
        
        // integration constant
        const double fac = intpoints.qwgt[iquad]*det;
        
        element_volume += fac;


        // get velocities (n+g,i) at integration point
        for (int isd=0;isd<NSD_;isd++)
        {
            velint[isd]=0.0;
            for (int inode=0;inode<numnode;inode++)
            {
                velint[isd] += funct[inode]*evelnp[isd+(NSD_*inode)];
            }
        }

        /*---------------- get history data (n,i) at integration point ---*/
        //expression for f3_veci(histvec,funct,evhist,iel);
        for (int isd=0;isd<NSD_;isd++)
        {
            histvec[isd]=0.0;
            for (int inode=0;inode<numnode;inode++)
            {
                histvec[isd] += funct[inode]*evhist[isd+(NSD_*inode)];
                //histvec[isd] += funct[inode]*evhist[isd][inode];
            }
        }

        // get velocity (np,i) derivatives at integration point
        for (int isd=0;isd<NSD_;isd++)
        {
            vderxy(0,isd)=0.0;
            vderxy(1,isd)=0.0;
            vderxy(2,isd)=0.0;
            for (int inode=0;inode<numnode;inode++)
            { 
                vderxy(0,isd) += derxy(isd,inode)*evelnp[0+(NSD_*inode)];
                vderxy(1,isd) += derxy(isd,inode)*evelnp[1+(NSD_*inode)];
                vderxy(2,isd) += derxy(isd,inode)*evelnp[2+(NSD_*inode)];
            }
        }

        // calculate 2nd velocity derivatives at integration point
        if (higher_order_ele)
        {
            for (int i=0;i<6;i++)
            {
                vderxy2(0,i)=0.0;
                vderxy2(1,i)=0.0;
                vderxy2(2,i)=0.0;
                for (int j=0;j<numnode;j++)
                {
                    vderxy2(0,i) += derxy2(i,j)*evelnp[0+(NSD_*j)];
                    vderxy2(1,i) += derxy2(i,j)*evelnp[1+(NSD_*j)];
                    vderxy2(2,i) += derxy2(i,j)*evelnp[2+(NSD_*j)];
                }
            }
        }



        // get grid velocity at integration point
        if(is_ale_)
            dserror("No ALE algorithms supported by Fluid3 element up to now.");
        else
        {
            gridvelint[0] = 0.0;
            gridvelint[1] = 0.0;
            gridvelint[2] = 0.0;
        }

        // get pressure gradients at integration point
        gradp[0] = gradp[1] = gradp[2] = 0.0;
        for (int inode=0; inode<numnode; inode++)
        {
            gradp[0] += derxy(0,inode) * eprenp[inode];
            gradp[1] += derxy(1,inode) * eprenp[inode];
            gradp[2] += derxy(2,inode) * eprenp[inode];            
        }
  
        // get pressure
        double press = 0.0;
        for (int inode=0;inode<numnode;inode++)
        {
            press += funct[inode]*eprenp[inode];
        }

        // get bodyforce in gausspoint
        for (int isd=0;isd<NSD_;isd++)
        {
            edeadng[isd] = 0.0;
            for (int i=0;i<numnode;i++)
            {
                edeadng[isd]+= bodyforce(isd,i)*funct[i];
            }
        }
      
        /*-------------- perform integration for entire matrix and rhs ---*/
        if(is_stationary==false)    
            f3_calmat(*sys_mat,*residual,velint,histvec,gridvelint,
                        press,vderxy,vderxy2,gradp,funct,tau,
                        derxy,derxy2,edeadng,fac,visc,numnode,params);
        else
            f3_calmat_stationary(*sys_mat,*residual,velint,histvec,gridvelint,
                        press,vderxy,vderxy2,gradp,funct,tau,
                        derxy,derxy2,edeadng,fac,visc,numnode,params);
     
    } // end of loop over integration points/

    //cout << "Element volume: " << element_volume << endl;

    return;
} // DRT::Elements::Fluid3::f3_sys_mat



/*----------------------------------------------------------------------*
 |  calculate Jacobian matrix                                           |
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix DRT::Elements::XFluid3::getJacobiMatrix(
                      const Epetra_SerialDenseMatrix& xyze,
                      const Epetra_SerialDenseMatrix& deriv,
                      const int                       iel) const
{
    
    Epetra_SerialDenseMatrix    xjm(NSD_,NSD_);
    
    // determine jacobian matrix at point r,s,t
    for (int isd=0; isd<NSD_; isd++)
    {
        for (int jsd=0; jsd<NSD_; jsd++)
        {
            double dum = 0.0;
            for (int inode=0; inode<iel; inode++)
            {
                dum += deriv(isd,inode)*xyze(jsd,inode);
            }
            xjm(isd,jsd) = dum;
        }
    }

    return xjm;
}


/*----------------------------------------------------------------------*
 |  calculate determinant                                               |
 *----------------------------------------------------------------------*/
double DRT::Elements::XFluid3::getDeterminante(const Epetra_SerialDenseMatrix&  xjm) const
{
    // determinant of jacobian matrix
    const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                       xjm(0,1)*xjm(1,2)*xjm(2,0)+
                       xjm(0,2)*xjm(1,0)*xjm(2,1)-
                       xjm(0,2)*xjm(1,1)*xjm(2,0)-
                       xjm(0,0)*xjm(1,2)*xjm(2,1)-
                       xjm(0,1)*xjm(1,0)*xjm(2,2);

    if(det < 0.0)
    {
        printf("\n");
        printf("GLOBAL ELEMENT NO.%i\n",Id());
        printf("NEGATIVE JACOBIAN DETERMINANT: %lf\n", det);
        dserror("Stopped not regulary!\n");
    }

    return det;
}






/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3::f3_getbodyforce(Epetra_SerialDenseMatrix& edeadng,
                                             const int                 iel    ,
                                             ParameterList&            params
)
{
    vector<DRT::Condition*> myneumcond;

    // check whether all nodes have a unique VolumeNeumann condition
    int count = 0;
    for(int nn=0;nn<iel;nn++)
    {
        Nodes()[nn]->GetCondition("VolumeNeumann",myneumcond);
        
        dsassert(myneumcond.size()<2, "more than one VolumeNeumann cond on one node");
        
        if (myneumcond.size()==1)
            count++;
    }

    if (count == iel)
    {
        // find out whether we will use a time curve
        bool usetime = true;
        const double time = params.get("total time",-1.0);
        if (time<0.0) usetime = false;

        const vector<int>* curve  = myneumcond[0]->Get<vector<int> >("curve");
        int curvenum = -1;

        // get the factor for the timecurve
        if (curve) curvenum = (*curve)[0];
        double curvefac = 1.0;
        if (curvenum>=0 && usetime)
            curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(time);

        // set this condition to the edeadng array
        for(int nn=0;nn<iel;nn++)
        {
            Nodes()[nn]->GetCondition("VolumeNeumann",myneumcond);

            // get values and switches from the condition
            const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
            const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

            for(int dim=0;dim<NSD_;dim++)
            {
                edeadng(dim,nn)=(*onoff)[dim]*(*val)[dim]*curvefac;
            }
        }
    }
    else
    {
        // we have no dead load
        for(int nn=0;nn<iel;nn++)
            {
            for(int dim=0;dim<NSD_;dim++)
                {
                edeadng(dim,nn)=0.0;
            }
        }
    }

    return;
}



/*----------------------------------------------------------------------*
 |  calculate global derivatives w.r.t. x,y,z at point r,s,t (private)genk05/02
 *----------------------------------------------------------------------*/
void DRT::Elements::XFluid3::f3_gder(Epetra_SerialDenseMatrix& derxy,
                                    const Epetra_SerialDenseMatrix& deriv,
                                    const Epetra_SerialDenseMatrix& xjm,
                                    const double& det,
                                    const int iel
                    ) const
{
  Epetra_SerialDenseMatrix  xji(NSD_,NSD_);   // inverse of jacobian matrix


 /*----------calculate global derivatives w.r.t. x,y,z at point r,s,t ---*/

    /*------------------------------------------------------- initialistion */
for(int k=0;k<iel;k++)
{
   derxy(0,k)=0.0;
   derxy(1,k)=0.0;
   derxy(2,k)=0.0;
} /* end of loop over k */

    /*------------------------------------------------- inverse of jacobian */
xji(0,0) = (  xjm(1,1)*xjm(2,2) - xjm(2,1)*xjm(1,2))/det;
xji(1,0) = (- xjm(1,0)*xjm(2,2) + xjm(2,0)*xjm(1,2))/det;
xji(2,0) = (  xjm(1,0)*xjm(2,1) - xjm(2,0)*xjm(1,1))/det;
xji(0,1) = (- xjm(0,1)*xjm(2,2) + xjm(2,1)*xjm(0,2))/det;
xji(1,1) = (  xjm(0,0)*xjm(2,2) - xjm(2,0)*xjm(0,2))/det;
xji(2,1) = (- xjm(0,0)*xjm(2,1) + xjm(2,0)*xjm(0,1))/det;
xji(0,2) = (  xjm(0,1)*xjm(1,2) - xjm(1,1)*xjm(0,2))/det;
xji(1,2) = (- xjm(0,0)*xjm(1,2) + xjm(1,0)*xjm(0,2))/det;
xji(2,2) = (  xjm(0,0)*xjm(1,1) - xjm(1,0)*xjm(0,1))/det;

    /*---------------------------------------- calculate global derivatives */
for (int k=0;k<iel;k++)
{
   derxy(0,k) +=   xji(0,0) * deriv(0,k) \
                  + xji(0,1) * deriv(1,k) \
                  + xji(0,2) * deriv(2,k) ;
   derxy(1,k) +=   xji(1,0) * deriv(0,k) \
                  + xji(1,1) * deriv(1,k) \
                  + xji(1,2) * deriv(2,k) ;
   derxy(2,k) +=   xji(2,0) * deriv(0,k) \
                  + xji(2,1) * deriv(1,k) \
                  + xji(2,2) * deriv(2,k) ;
} /* end of loop over k */

    /*----------------------------------------------------------------------*/

return;
} // end of DRT:Elements:Fluid3:f3_gder


void DRT::Elements::XFluid3::f3_gder2(const Epetra_SerialDenseMatrix& xyze,
                                      const Epetra_SerialDenseMatrix& xjm,
                                      const Epetra_SerialDenseMatrix& derxy,
                                      Epetra_SerialDenseMatrix& derxy2,
                                      const Epetra_SerialDenseMatrix& deriv2,
                                      const int iel
                                     ) const
{
double r0,r1,r2,r3,r4,r5;
//--------------------------------------------initialize and zero out everything
Epetra_SerialDenseMatrix bm(6,6);
Epetra_SerialDenseMatrix xder2(6,NSD_);

/*--------------------------- calculate elements of jacobian_bar matrix */
bm(0,0) = xjm(0,0)*xjm(0,0);
bm(1,0) = xjm(1,0)*xjm(1,0);
bm(2,0) = xjm(2,0)*xjm(2,0);
bm(3,0) = xjm(0,0)*xjm(1,0);
bm(4,0) = xjm(0,0)*xjm(2,0);
bm(5,0) = xjm(1,0)*xjm(2,0);

bm(0,1) = xjm(0,1)*xjm(0,1);
bm(1,1) = xjm(1,1)*xjm(1,1);
bm(2,1) = xjm(2,1)*xjm(2,1);
bm(3,1) = xjm(0,1)*xjm(1,1);
bm(4,1) = xjm(0,1)*xjm(2,1);
bm(5,1) = xjm(1,1)*xjm(2,1);

bm(0,2) = xjm(0,2)*xjm(0,2);
bm(1,2) = xjm(1,2)*xjm(1,2);
bm(2,2) = xjm(2,2)*xjm(2,2);
bm(3,2) = xjm(0,2)*xjm(1,2);
bm(4,2) = xjm(0,2)*xjm(2,2);
bm(5,2) = xjm(1,2)*xjm(2,2);

bm(0,3) = 2.0*xjm(0,0)*xjm(0,1);
bm(1,3) = 2.0*xjm(1,0)*xjm(1,1);
bm(2,3) = 2.0*xjm(2,0)*xjm(2,1);
bm(3,3) = xjm(0,0)*xjm(1,1)+xjm(1,0)*xjm(0,1);
bm(4,3) = xjm(0,0)*xjm(2,1)+xjm(2,0)*xjm(0,1);
bm(5,3) = xjm(1,0)*xjm(2,1)+xjm(2,0)*xjm(1,1);

bm(0,4) = 2.0*xjm(0,0)*xjm(0,2);
bm(1,4) = 2.0*xjm(1,0)*xjm(1,2);
bm(2,4) = 2.0*xjm(2,0)*xjm(2,2);
bm(3,4) = xjm(0,0)*xjm(1,2)+xjm(1,0)*xjm(0,2);
bm(4,4) = xjm(0,0)*xjm(2,2)+xjm(2,0)*xjm(0,2);
bm(5,4) = xjm(1,0)*xjm(2,2)+xjm(2,0)*xjm(1,2);

bm(0,5) = 2.0*xjm(0,1)*xjm(0,2);
bm(1,5) = 2.0*xjm(1,1)*xjm(1,2);
bm(2,5) = 2.0*xjm(2,1)*xjm(2,2);
bm(3,5) = xjm(0,1)*xjm(1,2)+xjm(1,1)*xjm(0,2);
bm(4,5) = xjm(0,1)*xjm(2,2)+xjm(2,1)*xjm(0,2);
bm(5,5) = xjm(1,1)*xjm(2,2)+xjm(2,1)*xjm(1,2);

/*-------------------------------------- inverse of jacobian_bar matrix */

LINALG::NonSymmetricInverse(bm,6);

// output for debug
// (for more details see the comments in definition of NonSymmetricInverse()
#if 0
for (int i = 0 ; i < 6; ++i)
{
for (int j = 0 ; j < 6; ++j)
{
       if (bm(i,j)!=0.0)
       printf("bm[%d][%d] %22.16e ",i,j,bm(i,j));
       else
        printf("bm[%d][%d] 0.000 ",i,j);
       printf("\n");
}
printf("\n");
}
#endif

/*----------------------------------------------------------- initialise*/
/*   already initialized by constructor of EpetraSerialDenseMeatrix
for (int i=0;i<NSD_;i++)
{
   for (int j=0;j<6;j++) xder2(j,i)=0.0;
}
*/

for (int i=0;i<iel;i++)
{
   for (int j=0;j<6;j++) derxy2(j,i)=0.0;
}

/*----------------------- determine 2nd derivatives of coord.-functions */
for (int i=0;i<iel;i++)
{
   xder2(0,0) += deriv2(0,i) * xyze(0,i);
   xder2(1,0) += deriv2(1,i) * xyze(0,i);
   xder2(2,0) += deriv2(2,i) * xyze(0,i);
   xder2(3,0) += deriv2(3,i) * xyze(0,i);
   xder2(4,0) += deriv2(4,i) * xyze(0,i);
   xder2(5,0) += deriv2(5,i) * xyze(0,i);

   xder2(0,1) += deriv2(0,i) * xyze(1,i);
   xder2(1,1) += deriv2(1,i) * xyze(1,i);
   xder2(2,1) += deriv2(2,i) * xyze(1,i);
   xder2(3,1) += deriv2(3,i) * xyze(1,i);
   xder2(4,1) += deriv2(4,i) * xyze(1,i);
   xder2(5,1) += deriv2(5,i) * xyze(1,i);

   xder2(0,2) += deriv2(0,i) * xyze(2,i);
   xder2(1,2) += deriv2(1,i) * xyze(2,i);
   xder2(2,2) += deriv2(2,i) * xyze(2,i);
   xder2(3,2) += deriv2(3,i) * xyze(2,i);
   xder2(4,2) += deriv2(4,i) * xyze(2,i);
   xder2(5,2) += deriv2(5,i) * xyze(2,i);
} /* end of loop over i */

/*--------------------------------- calculate second global derivatives */
for (int i=0;i<iel;i++)
{
   r0 = deriv2(0,i) - xder2(0,0)*derxy(0,i) - xder2(0,1)*derxy(1,i) \
                     - xder2(0,2)*derxy(2,i);
   r1 = deriv2(1,i) - xder2(1,0)*derxy(0,i) - xder2(1,1)*derxy(1,i) \
                     - xder2(1,2)*derxy(2,i);
   r2 = deriv2(2,i) - xder2(2,0)*derxy(0,i) - xder2(2,1)*derxy(1,i) \
                     - xder2(2,2)*derxy(2,i);
   r3 = deriv2(3,i) - xder2(3,0)*derxy(0,i) - xder2(3,1)*derxy(1,i) \
                     - xder2(3,2)*derxy(2,i);
   r4 = deriv2(4,i) - xder2(4,0)*derxy(0,i) - xder2(4,1)*derxy(1,i) \
                     - xder2(4,2)*derxy(2,i);
   r5 = deriv2(5,i) - xder2(5,0)*derxy(0,i) - xder2(5,1)*derxy(1,i) \
                     - xder2(5,2)*derxy(2,i);

   derxy2(0,i) += bm(0,0)*r0 + bm(0,1)*r1 + bm(0,2)*r2 \
                +  bm(0,3)*r3 + bm(0,4)*r4 + bm(0,5)*r5;
   derxy2(1,i) += bm(1,0)*r0 + bm(1,1)*r1 + bm(1,2)*r2 \
                +  bm(1,3)*r3 + bm(1,4)*r4 + bm(1,5)*r5;
   derxy2(2,i) += bm(2,0)*r0 + bm(2,1)*r1 + bm(2,2)*r2 \
                +  bm(2,3)*r3 + bm(2,4)*r4 + bm(2,5)*r5;
   derxy2(3,i) += bm(3,0)*r0 + bm(3,1)*r1 + bm(3,2)*r2 \
                +  bm(3,3)*r3 + bm(3,4)*r4 + bm(3,5)*r5;
   derxy2(4,i) += bm(4,0)*r0 + bm(4,1)*r1 + bm(4,2)*r2 \
                +  bm(4,3)*r3 + bm(4,4)*r4 + bm(4,5)*r5;
   derxy2(5,i) += bm(5,0)*r0 + bm(5,1)*r1 + bm(5,2)*r2 \
                +  bm(5,3)*r3 + bm(5,4)*r4 + bm(5,5)*r5;
} /* end of loop over i */

/*----------------------------------------------------------------------*/

return;
} // end of DRT:Elements:Fluid3:f3_gder2



/*----------------------------------------------------------------------*
 |  evaluate fluid coefficient matrix (private)              chfoe 04/04|
 *----------------------------------------------------------------------*/

/*
In this routine the Gauss point contributions to the elemental coefficient
matrix of a stabilised fluid3 element are calculated. The procedure is
based on the Rothe method of first integrating in time. Hence the
resulting terms include coefficients containing time integration variables
such as theta or delta t which are represented by 'timefac'.

The routine was completed to contain ALE-terms also.         chfoe 11/04

The stabilisation is based on the residuum:

R_M = u + timefac u * grad u - timefac * 2 nu div epsilon(u)
    + timefac grad p - rhsint

R_C = div u

The corresponding weighting operators are
L_M = v + timefac u_old * grad v + timefac v * grad u_old
    - timefac * 2 nu alpha div epsilon (v) + timefac beta grad q

L_C = div v

where alpha = -1
      beta  = -1
are sign regulating factors and rhsint differs for different time
These factores are worked in now and cannot be changed any more.

integration schemes:

One-step-Theta:
rhsint = u_old + Theta dt f + (1-Theta) acc_old

BDF2:

generalised alpha:


The stabilisation by means of the momentum residuum R_M is of the unusual
type:
   Galerkin parts MINUS sum over elements (stabilising parts)
The stabilisation by means of the continuity equation R_C is done in the
usual way:
   Galerkin parts PLUS sum over elements (stabilising parts)

The calculation proceeds as follows.
1) obtain single (linearised) operators of R_M, R_C, L_M and L_C
2) build Galerkin terms from them
3) build stabilising terms from them
4) build Galerkin and stabilising terms of RHS

NOTE: u_old represents the last iteration value. (The most recent one
      we've got!)

NOTE: Galerkin and stabilisation matrices are calculated within one
      routine.

NOTE: In order to increase the performance plenty of terms are concentrated
      and worked into each other. A lengthy version of the file is available
      from the author.


Notational remarks:

                   /              \
                  | u_x,x   u_x,y |
vderxy = grad u = |               |
                  | u_y,x   u_y,y |
                  \               /

           /                         \
          | u_x,xx   u_x,yy   u_x,xy |
vderxy2 = |                          |
          | u_y,xx   u_y,yy   u_y,xy |
          \                          /

for further comments see comment lines within code.

</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *eforce     DOUBLE        (o)   ele force vector
\param  *velint     DOUBLE        (i)   vel at INT point
\param  *histvec    DOUBLE        (i)   rhs at INT point
\param  *gridvint   DOUBLE        (i)   gridvel at INT point
\param **vderxy     DOUBLE        (i)   global vel derivatives
\param  *vderxy2    DOUBLE        (i)   2nd global vel derivatives
\param  *funct      DOUBLE        (i)   nat. shape funcs
\param **derxy      DOUBLE        (i)   global coord. deriv
\param **derxy2     DOUBLE        (i)   2nd global coord. deriv.
\param  *edeadng    DOUBLE        (i)   dead load at time n+1
\param   fac        DOUBLE        (i)   weighting factor
\param   visc       DOUBLE        (i)   fluid viscosity
\param   iel        INT           (i)   number of nodes of act. ele
\param  *hasext     INT           (i)   flag, if element has volume load
\param   isale      INT           (i)   flag, if ALE or EULER
\return void
------------------------------------------------------------------------*/

void DRT::Elements::XFluid3::f3_calmat( Epetra_SerialDenseMatrix& estif,
                Epetra_SerialDenseVector&  eforce,
                const vector<double>&            velint,
                const vector<double>&            histvec,
                const vector<double>&            gridvint,
                const double&                    press,
                const Epetra_SerialDenseMatrix&  vderxy,
                const Epetra_SerialDenseMatrix&  vderxy2,
                const vector<double>&            gradp,
                const Epetra_SerialDenseVector&  funct,
                const vector<double>&            tau,
                const Epetra_SerialDenseMatrix&  derxy,
                const Epetra_SerialDenseMatrix&  derxy2,
                const vector<double>&            edeadng,
                const double&                    fac,
                const double&                    visc,
                const int&                       iel,
                ParameterList&                   params
                )
{
    //DOUBLE  viscous[3][3][3*iel]; /* viscous term partially integrated */
    
    
    /*========================== initialisation ============================*/
    // One-step-Theta: timefac = theta*dt
    // BDF2:           timefac = 2/3 * dt
    const double timefac = params.get<double>("time constant for integration",-1.0);
    dsassert(timefac >= 0.0, "No time constant for integration supplied");
    
    // time step size
    //double dt = params.get<double>("delta time",-1.0);
    //  if (dt == -1.0) dserror("No dta supplied");
    
    // stabilisation parameter
    const double tau_M  = tau[0]*fac;
    const double tau_Mp = tau[1]*fac;
    const double tau_C  = tau[2]*fac;
    
    // integration factors and coefficients of single terms
    // double time2nue   = timefac * 2.0 * visc;
    const double timetauM   = timefac * tau_M;
    const double timetauMp  = timefac * tau_Mp;
    
    const double ttimetauM  = timefac * timetauM;
    const double ttimetauMp = timefac * timetauMp;
    const double timefacfac = timefac * fac;
    
    /*------------------------- evaluate rhs vector at integration point ---*/
    // no switch here at the moment w.r.t. is_ale
    vector<double>  rhsint(NSD_);  // total right hand side terms at int.-point
    rhsint[0] = histvec[0] + edeadng[0]*timefac;
    rhsint[1] = histvec[1] + edeadng[1]*timefac;
    rhsint[2] = histvec[2] + edeadng[2]*timefac;
    /*----------------- get numerical representation of single operators ---*/
    
    /* Convective term  u_old * grad u_old: */
    vector<double>  conv_old_temp(NSD_);
    for (int isd=0; isd<NSD_; isd++)
    {
        conv_old_temp[isd] = vderxy(isd,0) * velint[0] + vderxy(isd,1) * velint[1] + vderxy(isd,2) * velint[2];
    }
    const vector<double>  conv_old = conv_old_temp;
    
    /* new for incremental formulation: */
    /* Convective term  u_G_old * grad u_old: */
    vector<double>  conv_g_old_temp(NSD_);
    for (int isd=0; isd<NSD_; isd++)
    {
        conv_g_old_temp[isd] = (vderxy(isd,0) * gridvint[0] +
                                vderxy(isd,1) * gridvint[1] +
                                vderxy(isd,2) * gridvint[2]);
    }
    const vector<double>  conv_g_old = conv_g_old_temp;

    
    /* Viscous term  div epsilon(u_old) */
    vector<double>  visc_old(NSD_);
    visc_old[0] = vderxy2(0,0) + 0.5 * ( vderxy2(0,1) + vderxy2(1,3) + vderxy2(0,2) + vderxy2(2,4));
    visc_old[1] = vderxy2(1,1) + 0.5 * ( vderxy2(1,0) + vderxy2(0,3) + vderxy2(1,2) + vderxy2(2,5));
    visc_old[2] = vderxy2(2,2) + 0.5 * ( vderxy2(2,0) + vderxy2(0,4) + vderxy2(2,1) + vderxy2(1,5));
    
    
    vector<double>  conv_c(iel);        /* linearisation of convect, convective part */
    vector<double>  conv_g(iel);        /* linearisation of convect, grid part */
    Epetra_SerialDenseMatrix  ugradv(iel,NSD_*iel);    /* linearisation of u * grad v   */
    Epetra_SerialDenseMatrix  vconv_r(NSD_,iel);
    Epetra_SerialDenseMatrix  conv_r(NSD_,NSD_*iel);  /* linearisation of convect, reactive part */
    Epetra_SerialDenseMatrix  viscs2(NSD_,NSD_*iel);      /* viscous term incluiding 2nd derivatives */
    vector<double>  div(NSD_*iel);             /* divergence of u or v              */
    for (int inode=0; inode<iel; inode++) /* loop over nodes of element */
    {
       /* Reactive term  u:  funct */
       /* linearise convective term */
    
       /*--- convective part u_old * grad (funct) --------------------------*/
       /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
          with  N .. form function matrix                                   */
       conv_c[inode] = derxy(0,inode) * velint[0] 
                     + derxy(1,inode) * velint[1]
                     + derxy(2,inode) * velint[2];
    
       /*--- convective grid part u_G * grad (funct) -----------------------*/
       /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
       if(is_ale_)
       {
         dserror("No ALE supported by Fluid3 at the moment.");
          //    conv_g[i] = - derxy(0,i) * gridvint[0] - derxy(1,i) * gridvint[1]
          //           - derxy(2,i) * gridvint[2];
       }
       else
       {
         conv_g[inode] = 0.0;
       }
    
       /*--- reactive part funct * grad (u_old) ----------------------------*/
       /* /                                     \
          |  u_old_x,x   u_old_x,y   u_old x,z  |
          |                                     |
          |  u_old_y,x   u_old_y,y   u_old_y,z  | * N
          |                                     |
          |  u_old_z,x   u_old_z,y   u_old_z,z  |
          \                                     /
          with  N .. form function matrix                                   */
    
       conv_r(0,NSD_*inode)   = vderxy(0,0)*funct[inode];
       conv_r(0,NSD_*inode+1) = vderxy(0,1)*funct[inode];
       conv_r(0,NSD_*inode+2) = vderxy(0,2)*funct[inode];
       conv_r(1,NSD_*inode)   = vderxy(1,0)*funct[inode];
       conv_r(1,NSD_*inode+1) = vderxy(1,1)*funct[inode];
       conv_r(1,NSD_*inode+2) = vderxy(1,2)*funct[inode];
       conv_r(2,NSD_*inode)   = vderxy(2,0)*funct[inode];
       conv_r(2,NSD_*inode+1) = vderxy(2,1)*funct[inode];
       conv_r(2,NSD_*inode+2) = vderxy(2,2)*funct[inode];
    
       vconv_r(0,inode) = conv_r(0,NSD_*inode)*velint[0] + conv_r(0,NSD_*inode+1)*velint[1] + conv_r(0,NSD_*inode+2)*velint[2];
       vconv_r(1,inode) = conv_r(1,NSD_*inode)*velint[0] + conv_r(1,NSD_*inode+1)*velint[1] + conv_r(1,NSD_*inode+2)*velint[2];
       vconv_r(2,inode) = conv_r(2,NSD_*inode)*velint[0] + conv_r(2,NSD_*inode+1)*velint[1] + conv_r(2,NSD_*inode+2)*velint[2];
    
       /*--- viscous term  - grad * epsilon(u): ----------------------------*/
       /*   /                                                \
            |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
          1 |                                                |
        - - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
          2 |                                                |
            |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
            \                                                /
    
        with N_x .. x-line of N
             N_y .. y-line of N                                             */
    
       viscs2(0,NSD_*inode)   = - 0.5 * (2.0 * derxy2(0,inode) + derxy2(1,inode) + derxy2(2,inode));
       viscs2(0,NSD_*inode+1) = - 0.5 *  derxy2(3,inode);
       viscs2(0,NSD_*inode+2) = - 0.5 *  derxy2(4,inode);
       viscs2(1,NSD_*inode)   = - 0.5 *  derxy2(3,inode);
       viscs2(1,NSD_*inode+1) = - 0.5 * (derxy2(0,inode) + 2.0 * derxy2(1,inode) + derxy2(2,inode));
       viscs2(1,NSD_*inode+2) = - 0.5 *  derxy2(5,inode);
       viscs2(2,NSD_*inode)   = - 0.5 *  derxy2(4,inode);
       viscs2(2,NSD_*inode+1) = - 0.5 *  derxy2(5,inode);
       viscs2(2,NSD_*inode+2) = - 0.5 * (derxy2(0,inode) + derxy2(1,inode) + 2.0 * derxy2(2,inode));
    
       /*--- viscous term (after integr. by parts) -------------------------*/
       /*   /                                             \
            |  2 N_x,x    N_x,y + N_y,x    N_x,z + N_z,x  |
          1 |                                             |
          - |  N_y,x + N_x,y   2 N_y,y     N_y,z + N_z,y  |
          2 |                                             |
            |  N_z,x + N_x,z   N_z,y + N_y,z    2 N_z,z   |
            \                                             /
       with N_x .. x-line of N
            N_y .. y-line of N
            N_z .. z-line of N                                              */


    /* not needed for incremental solver routine     g.bau 03/07
       viscous[0][0][3*i]   = derxy(0,i);
       viscous[0][0][3*i+1] = 0.0;
       viscous[0][0][3*i+2] = 0.0;                // 1st index:
       viscous[0][1][3*i]   = 0.5 * derxy(1,i);  //   line of epsilon
       viscous[0][1][3*i+1] = 0.5 * derxy(0,i);  // 2nd index:
       viscous[0][1][3*i+2] = 0.0;                //   column of epsilon
       viscous[0][2][3*i]   = 0.5 * derxy(2,i);  // 3rd index:
       viscous[0][2][3*i+1] = 0.0;                //   elemental vel dof
       viscous[0][2][3*i+2] = 0.5 * derxy(0,i);
       viscous[1][0][3*i]   = 0.5 * derxy(1,i);
       viscous[1][0][3*i+1] = 0.5 * derxy(0,i);
       viscous[1][0][3*i+2] = 0.0;
       viscous[1][1][3*i]   = 0.0;
       viscous[1][1][3*i+1] = derxy(1,i);
       viscous[1][1][3*i+2] = 0.0;
       viscous[1][2][3*i]   = 0.0;
       viscous[1][2][3*i+1] = 0.5 * derxy(2,i);
       viscous[1][2][3*i+2] = 0.5 * derxy(1,i);
       viscous[2][0][3*i]   = 0.5 * derxy(2,i);
       viscous[2][0][3*i+1] = 0.0;
       viscous[2][0][3*i+2] = 0.5 * derxy(0,i);
       viscous[2][1][3*i]   = 0.0;
       viscous[2][1][3*i+1] = 0.5 * derxy(2,i);
       viscous[2][1][3*i+2] = 0.5 * derxy(1,i);
       viscous[2][2][3*i]   = 0.0;
       viscous[2][2][3*i+1] = 0.0;
       viscous[2][2][3*i+2] = derxy(2,i);
    */
    
       /* pressure gradient term derxy, funct without or with integration   *
        * by parts, respectively                                            */
    
       /*--- divergence u term ---------------------------------------------*/
       div[NSD_*inode]   = derxy(0,inode);
       div[NSD_*inode+1] = derxy(1,inode);
       div[NSD_*inode+2] = derxy(2,inode);
    
       /*--- ugradv-Term ---------------------------------------------------*/
       /*
         /                                                          \
         |  N1*N1,x  N1*N1,y  N2*N1,x  N2*N1,y  N3*N1,x ...       . |
         |                                                          |
         |  N1*N2,x  N1*N2,y  N2*N2,x  N2*N2,y  N3*N2,x ...       . |
         |                                                          |
         |  N1*N3,x  N1*N3,y  N2*N3,x  N2*N3,y  N3*N3,x ...       . |
         |                                           .              |
         |  . . .                                        .          |
         |                                                  Ni*Ni,y |
         \                                                          /       */
       /* remark: vgradu = ugradv^T */
       for (int j=0; j<iel; j++)
       {
          ugradv(inode,NSD_*j)   = derxy(0,inode) * funct[j];
          ugradv(inode,NSD_*j+1) = derxy(1,inode) * funct[j];
          ugradv(inode,NSD_*j+2) = derxy(2,inode) * funct[j];
       }
    
    } // end of loop over nodes of element

    /*--------------------------------- now build single stiffness terms ---*/
    
    #define estif_(i,j)    estif(i,j)
    #define eforce_(i)     eforce[i]
    #define funct_(i)      funct[i]
    #define vderxyz_(i,j)  vderxy(i,j)
    #define conv_c_(j)     conv_c[j]
    #define conv_g_(j)     conv_g[j]
    #define conv_r_(i,j,k) conv_r(i,3*(k)+j)
    #define vconv_r_(i,j)  vconv_r(i,j)
    #define conv_old_(j)   conv_old[j]
    #define conv_g_old_(j) conv_g_old[j]
    #define derxyz_(i,j)   derxy(i,j)
    #define gridvint_(j)   gridvint[j]
    #define velint_(j)     velint[j]
    #define viscs2_(i,j,k) viscs2(i,3*(k)+j)
    #define visc_old_(i)   visc_old[i]
    #define rhsint_(i)     rhsint[i]
    #define gradp_(j)      gradp[j]
    #define nu_            visc
    #define thsl           timefac
    
      /* This code is generated using MuPAD. Ask me for the MuPAD. u.kue */
    
      /* We keep two versions: with and without ale. The laster one is a
       * little faster. (more than 10%) */
    
      //if (!is_ale_)
      {
        #include "fluid3_xfem_stiff.cpp"
        #include "fluid3_xfem_rhs_incr.cpp"
      }
      /*
      else
      {
      dserror("No ALE support in Fluid3.");
      // #include "f3_stiff_ale.c"
      }
      */
    
    #undef estif_
    #undef eforce_
    #undef conv_c_
    #undef conv_g_
    #undef conv_r_
    #undef vconv_r_
    #undef conv_old_
    #undef conv_g_old_
    #undef derxyz_
    #undef gridvint_
    #undef velint_
    #undef viscs2_
    #undef gradp_
    #undef funct_
    #undef vderxyz_
    #undef visc_old_
    #undef rhsint_
    #undef nu_
    #undef thsl
    
    return;
} // end of DRT:Elements:Fluid3:f3_calmat


void DRT::Elements::XFluid3::f3_calmat_stationary( Epetra_SerialDenseMatrix& estif,
                Epetra_SerialDenseVector&  eforce,
                const vector<double>&            velint,
                const vector<double>&            histvec,
                const vector<double>&            gridvint,
                const double&                    press,
                const Epetra_SerialDenseMatrix&  vderxy,
                const Epetra_SerialDenseMatrix&  vderxy2,
                const vector<double>&            gradp,
                const Epetra_SerialDenseVector&  funct,
                const vector<double>&            tau,
                const Epetra_SerialDenseMatrix&  derxy,
                const Epetra_SerialDenseMatrix&  derxy2,
                const vector<double>&            edeadng,
                const double&                    fac,
                const double&                    visc,
                const int&                       iel,
                ParameterList&                   params
                )
{

/*========================= further variables =========================*/

Epetra_SerialDenseMatrix  viscs2(NSD_,NSD_*iel);      /* viscous term incluiding 2nd derivatives */
vector<double>  conv_c(iel);        /* linearisation of convect, convective part */
vector<double>  conv_g(iel);        /* linearisation of convect, grid part */
Epetra_SerialDenseMatrix  conv_r(NSD_,NSD_*iel);  /* linearisation of convect, reactive part */
vector<double>  div(NSD_*iel);             /* divergence of u or v              */
Epetra_SerialDenseMatrix  ugradv(iel,NSD_*iel);    /* linearisation of u * grad v   */
vector<double>  conv_old(NSD_);        /* convective term evalaluated with old velocities */
vector<double>  conv_g_old(NSD_);
vector<double>  visc_old(NSD_);        /* viscous term evaluated with old velocities      */
vector<double>  rhsint(NSD_);          /* total right hand side terms at int.-point       */
Epetra_SerialDenseMatrix  vconv_r(NSD_,iel);

/*========================== initialisation ============================*/
// One-step-Theta: timefac = theta*dt
// BDF2:           timefac = 2/3 * dt
  const double timefac = params.get<double>("time constant for integration",-1.0);
  if (timefac < 0.0) dserror("No time constant for integration supplied");

// stabilisation parameter
const double tau_M  = tau[0]*fac;
const double tau_Mp = tau[1]*fac;
const double tau_C  = tau[2]*fac;

/*------------------------- evaluate rhs vector at integration point ---*/
// no switch here at the moment w.r.t. is_ale
    rhsint[0] = histvec[0] + edeadng[0]*timefac;
    rhsint[1] = histvec[1] + edeadng[1]*timefac;
    rhsint[2] = histvec[2] + edeadng[2]*timefac;
/*----------------- get numerical representation of single operators ---*/

/* Convective term  u_old * grad u_old: */
conv_old[0] = vderxy(0,0) * velint[0] + vderxy(0,1) * velint[1]
            + vderxy(0,2) * velint[2];
conv_old[1] = vderxy(1,0) * velint[0] + vderxy(1,1) * velint[1]
            + vderxy(1,2) * velint[2];
conv_old[2] = vderxy(2,0) * velint[0] + vderxy(2,1) * velint[1]
            + vderxy(2,2) * velint[2];

/* new for incremental formulation: */
/* Convective term  u_G_old * grad u_old: */
conv_g_old[0] = (vderxy(0,0) * gridvint[0] +
         vderxy(0,1) * gridvint[1] +
         vderxy(0,2) * gridvint[2]);
conv_g_old[1] = (vderxy(1,0) * gridvint[0] +
         vderxy(1,1) * gridvint[1] +
         vderxy(1,2) * gridvint[2]);
conv_g_old[2] = (vderxy(2,0) * gridvint[0] +
         vderxy(2,1) * gridvint[1] +
         vderxy(2,2) * gridvint[2]);

/* Viscous term  div epsilon(u_old) */
visc_old[0] = vderxy2(0,0) + 0.5 * ( vderxy2(0,1) + vderxy2(1,3)
                                    + vderxy2(0,2) + vderxy2(2,4));
visc_old[1] = vderxy2(1,1) + 0.5 * ( vderxy2(1,0) + vderxy2(0,3)
                                    + vderxy2(1,2) + vderxy2(2,5));
visc_old[2] = vderxy2(2,2) + 0.5 * ( vderxy2(2,0) + vderxy2(0,4)
                                    + vderxy2(2,1) + vderxy2(1,5));

for (int i=0; i<iel; i++) /* loop over nodes of element */
{
   /* Reactive term  u:  funct */
   /* linearise convective term */

   /*--- convective part u_old * grad (funct) --------------------------*/
   /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
      with  N .. form function matrix                                   */
   conv_c[i] = derxy(0,i) * velint[0] + derxy(1,i) * velint[1]
             + derxy(2,i) * velint[2];

   /*--- convective grid part u_G * grad (funct) -----------------------*/
   /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
   if(is_ale_)
   {
     dserror("No ALE supported by Fluid3 at the moment.");
      //    conv_g[i] = - derxy(0,i) * gridvint[0] - derxy(1,i) * gridvint[1]
      //           - derxy(2,i) * gridvint[2];
   }
   else
   {
     conv_g[i] = 0.0;
   }

   /*--- reactive part funct * grad (u_old) ----------------------------*/
   /* /                                     \
      |  u_old_x,x   u_old_x,y   u_old x,z  |
      |                                     |
      |  u_old_y,x   u_old_y,y   u_old_y,z  | * N
      |                                     |
      |  u_old_z,x   u_old_z,y   u_old_z,z  |
      \                                     /
      with  N .. form function matrix                                   */

   conv_r(0,NSD_*i)   = vderxy(0,0)*funct[i];
   conv_r(0,NSD_*i+1) = vderxy(0,1)*funct[i];
   conv_r(0,NSD_*i+2) = vderxy(0,2)*funct[i];
   conv_r(1,NSD_*i)   = vderxy(1,0)*funct[i];
   conv_r(1,NSD_*i+1) = vderxy(1,1)*funct[i];
   conv_r(1,NSD_*i+2) = vderxy(1,2)*funct[i];
   conv_r(2,NSD_*i)   = vderxy(2,0)*funct[i];
   conv_r(2,NSD_*i+1) = vderxy(2,1)*funct[i];
   conv_r(2,NSD_*i+2) = vderxy(2,2)*funct[i];

   vconv_r(0,i) = conv_r(0,NSD_*i)*velint[0] + conv_r(0,NSD_*i+1)*velint[1] + conv_r(0,NSD_*i+2)*velint[2];
   vconv_r(1,i) = conv_r(1,NSD_*i)*velint[0] + conv_r(1,NSD_*i+1)*velint[1] + conv_r(1,NSD_*i+2)*velint[2];
   vconv_r(2,i) = conv_r(2,NSD_*i)*velint[0] + conv_r(2,NSD_*i+1)*velint[1] + conv_r(2,NSD_*i+2)*velint[2];

   /*--- viscous term  - grad * epsilon(u): ----------------------------*/
   /*   /                                                \
        |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
      1 |                                                |
    - - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
      2 |                                                |
        |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
        \                                                /

    with N_x .. x-line of N
         N_y .. y-line of N                                             */

   viscs2(0,NSD_*i)   = - 0.5 * (2.0 * derxy2(0,i) + derxy2(1,i) + derxy2(2,i));
   viscs2(0,NSD_*i+1) = - 0.5 *  derxy2(3,i);
   viscs2(0,NSD_*i+2) = - 0.5 *  derxy2(4,i);
   viscs2(1,NSD_*i)   = - 0.5 *  derxy2(3,i);
   viscs2(1,NSD_*i+1) = - 0.5 * (derxy2(0,i) + 2.0 * derxy2(1,i) + derxy2(2,i));
   viscs2(1,NSD_*i+2) = - 0.5 *  derxy2(5,i);
   viscs2(2,NSD_*i)   = - 0.5 *  derxy2(4,i);
   viscs2(2,NSD_*i+1) = - 0.5 *  derxy2(5,i);
   viscs2(2,NSD_*i+2) = - 0.5 * (derxy2(0,i) + derxy2(1,i) + 2.0 * derxy2(2,i));

   /*--- viscous term (after integr. by parts) -------------------------*/
   /*   /                                             \
        |  2 N_x,x    N_x,y + N_y,x    N_x,z + N_z,x  |
      1 |                                             |
      - |  N_y,x + N_x,y   2 N_y,y     N_y,z + N_z,y  |
      2 |                                             |
        |  N_z,x + N_x,z   N_z,y + N_y,z    2 N_z,z   |
        \                                             /
   with N_x .. x-line of N
        N_y .. y-line of N
        N_z .. z-line of N                                              */


   /* pressure gradient term derxy, funct without or with integration   *
    * by parts, respectively                                            */

   /*--- divergence u term ---------------------------------------------*/
   div[NSD_*i]   = derxy(0,i);
   div[NSD_*i+1] = derxy(1,i);
   div[NSD_*i+2] = derxy(2,i);

   /*--- ugradv-Term ---------------------------------------------------*/
   /*
     /                                                          \
     |  N1*N1,x  N1*N1,y  N2*N1,x  N2*N1,y  N3*N1,x ...       . |
     |                                                          |
     |  N1*N2,x  N1*N2,y  N2*N2,x  N2*N2,y  N3*N2,x ...       . |
     |                                                          |
     |  N1*N3,x  N1*N3,y  N2*N3,x  N2*N3,y  N3*N3,x ...       . |
     |                                           .              |
     |  . . .                                        .          |
     |                                                  Ni*Ni,y |
     \                                                          /       */
   /* remark: vgradu = ugradv^T */
   for (int j=0; j<iel; j++)
   {
      ugradv(i,NSD_*j)   = derxy(0,i) * funct[j];
      ugradv(i,NSD_*j+1) = derxy(1,i) * funct[j];
      ugradv(i,NSD_*j+2) = derxy(2,i) * funct[j];
   }

} // end of loop over nodes of element

/*--------------------------------- now build single stiffness terms ---*/

#define estif_(i,j)    estif(i,j)
#define eforce_(i)     eforce[i]
#define funct_(i)      funct[i]
#define vderxyz_(i,j)  vderxy(i,j)
#define conv_c_(j)     conv_c[j]
#define conv_g_(j)     conv_g[j]
#define conv_r_(i,j,k) conv_r(i,3*(k)+j)
#define vconv_r_(i,j)  vconv_r(i,j)
#define conv_old_(j)   conv_old[j]
#define conv_g_old_(j) conv_g_old[j]
#define derxyz_(i,j)   derxy(i,j)
#define gridvint_(j)   gridvint[j]
#define velint_(j)     velint[j]
#define viscs2_(i,j,k) viscs2(i,3*(k)+j)
#define visc_old_(i)   visc_old[i]
#define rhsint_(i)     rhsint[i]
#define gradp_(j)      gradp[j]
#define nu_            visc
#define thsl           timefac

  /* This code is generated using MuPAD. Ask me for the MuPAD. u.kue */

  {
    #include "fluid3_xfem_stiff_stationary.cpp"
    #include "fluid3_xfem_rhs_incr_stationary.cpp"
  }

#undef estif_
#undef eforce_
#undef conv_c_
#undef conv_g_
#undef conv_r_
#undef vconv_r_
#undef conv_old_
#undef conv_g_old_
#undef derxyz_
#undef gridvint_
#undef velint_
#undef viscs2_
#undef gradp_
#undef funct_
#undef vderxyz_
#undef visc_old_
#undef rhsint_
#undef nu_
#undef thsl

return;
} // end of DRT:Elements:Fluid3:f3_calmat_stationary


// get optimal gaussrule for discretization type
GaussRule3D DRT::Elements::XFluid3::getOptimalGaussrule(const DiscretizationType& distype)
{
    GaussRule3D rule;
    switch (distype)
    {
    case hex8:
        rule = intrule_hex_8point;
        break;
    case hex20: case hex27:
        rule = intrule_hex_27point;
        break;
    case tet4:
        rule = intrule_tet_4point;
        break;
    case tet10:
        rule = intrule_tet_10point;
        break;
    default: 
        dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}

/*---------------------------------------------------------------------*
 |  calculate error for beltrami test problem (private)     gammi 04/07|
 *---------------------------------------------------------------------*/
void DRT::Elements::XFluid3::f3_int_beltrami_err(
  const vector<double>&           evelnp,
  const vector<double>&           eprenp,
  const struct _MATERIAL*         material,
  ParameterList&        params
  )
{
    /*-------------------------- add element error to "integrated" error */
    double velerr = params.get<double>("L2 integrated velocity error");
    double preerr = params.get<double>("L2 integrated pressure error");

    /*------------------------------------------------- set element data */
    const int iel = NumNode();
    const DiscretizationType distype = this->Shape();

    Epetra_SerialDenseVector  funct(iel);
    Epetra_SerialDenseMatrix  xjm(NSD_,NSD_);
    Epetra_SerialDenseMatrix  deriv(NSD_,iel);

    // get node coordinates of element
    Epetra_SerialDenseMatrix xyze(NSD_,iel);
    for(int inode=0;inode<iel;inode++)
    {
        xyze(0,inode)=Nodes()[inode]->X()[0];
        xyze(1,inode)=Nodes()[inode]->X()[1];
        xyze(2,inode)=Nodes()[inode]->X()[2];
    }

    //------------------------------ set constants for analytical solution
    const double t = params.get("total time",-1.0);
    if (t<0)
    {
        dserror("beltrami: no total time for error calculation");
    }

    const double a      = PI/4.0;
    const double d      = PI/2.0;

    /* get viscosity ---*/
    const double  visc = material->m.fluid->viscosity;


    const IntegrationPoints3D  intpoint = getIntegrationPoints3D(gaussrule_);

    vector<double> velint  (NSD_);
    vector<double> xint    (NSD_);

    
    vector<double> u       (NSD_);

    vector<double> deltavel(NSD_);

    // start loop over integration points
    for (int iquad=0;iquad<intpoint.nquad;iquad++)
    {
        // get values of  shape functions and their derivatives
        const double e1   = intpoint.qxg[iquad][0];
        const double e2   = intpoint.qxg[iquad][1];
        const double e3   = intpoint.qxg[iquad][2];
        shape_function_3D(funct ,e1,e2,e3,distype);
        shape_function_3D_deriv1(deriv ,e1,e2,e3,distype);

        // compute Jacobian matrix
        const Epetra_SerialDenseMatrix xjm = getJacobiMatrix(xyze,deriv,iel);
        const double det = getDeterminante(xjm);
        const double fac = intpoint.qwgt[iquad]*det;

        // get velocity at integration point
        for (int isd=0; isd<NSD_; isd++)
        {
            velint[isd] = 0.0;
            for (int inode=0;inode<iel;inode++)
            {
                velint[isd] += funct[inode]*evelnp[isd+(NSD_*inode)];
            }
        }

        // get pressure sol at integration point
        double preint = 0.0;
        for (int inode=0;inode<iel;inode++)
        {
            preint += funct[inode]*eprenp[inode];
        }

        // get velocity sol at integration point
        for (int isd=0; isd<NSD_; isd++)
        {
            xint[isd] = 0.0;
            for (int inode=0;inode<iel;inode++)
            {
                xint[isd] += funct[inode]*xyze(isd,inode);
            }
        }


        // compute analytical pressure
        const double p = -a*a/2.0 *
              ( exp(2.0*a*xint[0])
              + exp(2.0*a*xint[1])
              + exp(2.0*a*xint[2])
              + 2.0 * sin(a*xint[0] + d*xint[1]) * cos(a*xint[2] + d*xint[0]) * exp(a*(xint[1]+xint[2]))
              + 2.0 * sin(a*xint[1] + d*xint[2]) * cos(a*xint[0] + d*xint[1]) * exp(a*(xint[2]+xint[0]))
              + 2.0 * sin(a*xint[2] + d*xint[0]) * cos(a*xint[1] + d*xint[2]) * exp(a*(xint[0]+xint[1]))
              )* exp(-2.0*visc*d*d*t);

        // compute analytical velocities
        u[0] = -a * ( exp(a*xint[0]) * sin(a*xint[1] + d*xint[2]) +
                      exp(a*xint[2]) * cos(a*xint[0] + d*xint[1]) ) * exp(-visc*d*d*t);
        u[1] = -a * ( exp(a*xint[1]) * sin(a*xint[2] + d*xint[0]) +
                      exp(a*xint[0]) * cos(a*xint[1] + d*xint[2]) ) * exp(-visc*d*d*t);
        u[2] = -a * ( exp(a*xint[2]) * sin(a*xint[0] + d*xint[1]) +
                      exp(a*xint[1]) * cos(a*xint[2] + d*xint[0]) ) * exp(-visc*d*d*t);

        // compute difference between analytical solution and numerical solution
        const double deltap = preint - p;

        for (int isd=0;isd<NSD_;isd++)
        {
            deltavel[isd]=velint[isd]-u[isd];
        }

        // add square to L2 error
        for (int isd=0;isd<NSD_;isd++)
        {
            velerr += deltavel[isd]*deltavel[isd]*fac;
        }
        preerr += deltap*deltap*fac;

    } // end of loop over integration points


    // we use the parameterlist as a container to transport the calculated
    // errors from the elements to the dynamic routine

    params.set<double>("L2 integrated velocity error",velerr);
    params.set<double>("L2 integrated pressure error",preerr);

    return;
}



vector<double> DRT::Elements::XFluid3::f3_caltau(
    const Epetra_SerialDenseMatrix&         xyze,
    const vector<double>&                   evelnp,
    const DRT::Element::DiscretizationType  distype,
    const double                            visc,
    const int                               numnode,
    const double                            timefac,
    const bool                              is_stationary
    ) const
{
    // use one point gauss rule to calculate tau at element center
    GaussRule3D integrationrule_stabili;
    switch(distype)
    {
    case hex8: case hex20: case hex27:
        integrationrule_stabili = intrule_hex_1point;
        break;
    case tet4: case tet10:
        integrationrule_stabili = intrule_tet_1point;
        break;
    default: 
        dserror("invalid discretization type for fluid3");
    }

    // gaussian points
    const IntegrationPoints3D  intpoints = getIntegrationPoints3D(integrationrule_stabili);

    // shape functions and derivs at element center
    const double e1    = intpoints.qxg[0][0];
    const double e2    = intpoints.qxg[0][1];
    const double e3    = intpoints.qxg[0][2];
    const double wquad = intpoints.qwgt[0];
    
    Epetra_SerialDenseVector    funct(numnode);
    Epetra_SerialDenseMatrix    deriv(NSD_, numnode);
    shape_function_3D(funct,e1,e2,e3,distype);
    shape_function_3D_deriv1(deriv,e1,e2,e3,distype);
    
    // get element type constant for tau
    double mk=0.0;
    switch(distype)
    {
    case tet4: case hex8:
        mk = 0.333333333333333333333;
        break;
    case hex20: case hex27: case tet10:
        mk = 0.083333333333333333333;
        break;
    default: 
        dserror("type unknown!\n");
    }
    
    // get velocities at element center
    vector<double>  velint(NSD_);
    for (int isd=0;isd<NSD_;isd++)
    {
        velint[isd]=0.0;
        for (int inode=0;inode<numnode;inode++)
        {
            velint[isd] += funct[inode]*evelnp[isd+(NSD_*inode)];
        }
    }

    // get Jacobian matrix and determinant
    const Epetra_SerialDenseMatrix xjm = getJacobiMatrix(xyze,deriv,numnode);
    const double det = getDeterminante(xjm);
    const double vol = wquad*det;

    // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
    const double hk = pow((SIX*vol/PI),(1.0/3.0))/sqrt(3.0);

    // get derivatives
    Epetra_SerialDenseMatrix    derxy(NSD_, numnode);
    f3_gder(derxy,deriv,xjm,det,numnode);

    // get velocity norm
    const double vel_norm=sqrt( velint[0]*velint[0]
                              + velint[1]*velint[1]
                              + velint[2]*velint[2]);
    
    // normed velocity at element centre
    vector<double>  velino(NSD_);     
    if(vel_norm>=EPS6)
    {
        velino[0] = velint[0]/vel_norm;
        velino[1] = velint[1]/vel_norm;
        velino[2] = velint[2]/vel_norm;
    }
    else
    {
        velino[0] = 1.0;
        velino[1] = 0.0;
        velino[2] = 0.0;
    }
    
    // get streamlength
    double val = 0.0;
    for (int inode=0;inode<numnode;inode++)
    {
        val += abs(velino[0]*derxy(0,inode) 
                  +velino[1]*derxy(1,inode) 
                  +velino[2]*derxy(2,inode));
    }
    const double strle = 2.0/val;

    // calculate tau
    vector<double>  tau(3); // stab parameters
    if (is_stationary == false)
    {// stabilization parameters for instationary case (default)

        /*----------------------------------------------------- compute tau_Mu ---*/
        /* stability parameter definition according to

                  Barrenechea, G.R. and Valentin, F.: An unusual stabilized finite
                  element method for a generalized Stokes problem. Numerische
                  Mathematik, Vol. 92, pp. 652-677, 2002.
                  http://www.lncc.br/~valentin/publication.htm
        and:
                  Franca, L.P. and Valentin, F.: On an Improved Unusual Stabilized
                  Finite Element Method for the Advective-Reactive-Diffusive
                  Equation. Computer Methods in Applied Mechanics and Enginnering,
                  Vol. 190, pp. 1785-1800, 2000.
                  http://www.lncc.br/~valentin/publication.htm                   */


        const double re1 =/* 2.0*/ 4.0 * timefac * visc / (mk * DSQR(strle)); /* viscous : reactive forces */
        const double re2 = mk * vel_norm * strle / /* *1.0 */(2.0 * visc);    /* convective : viscous forces */

        const double xi1 = DMAX(re1,1.0);
        const double xi2 = DMAX(re2,1.0);

        tau[0] = DSQR(strle) / (DSQR(strle)*xi1+(/* 2.0*/ 4.0 * timefac*visc/mk)*xi2);

        // compute tau_Mp
        //    stability parameter definition according to Franca and Valentin (2000)
        //                                       and Barrenechea and Valentin (2002)
        const double re_viscous = /* 2.0*/ 4.0 * timefac * visc / (mk * DSQR(hk)); /* viscous : reactive forces */
        const double re_convect = mk * vel_norm * hk / /* *1.0 */(2.0 * visc);     /* convective : viscous forces */

        const double xi_viscous = DMAX(re_viscous,1.0);
        const double xi_convect = DMAX(re_convect,1.0);

        /*
                  xi1,xi2 ^
                          |      /
                          |     /
                          |    /
                        1 +---+
                          |
                          |
                          |
                          +--------------> re1,re2
                              1
        */
        tau[1] = DSQR(hk) / (DSQR(hk) * xi_viscous + (/* 2.0*/ 4.0 * timefac * visc/mk) * xi_convect);
    
        /*------------------------------------------------------ compute tau_C ---*/
        /*-- stability parameter definition according to Codina (2002), CMAME 191
         *
         * Analysis of a stabilized finite element approximation of the transient
         * convection-diffusion-reaction equation using orthogonal subscales.
         * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
         *
         * */
        //tau[2] = sqrt(DSQR(visc)+DSQR(0.5*vel_norm*hk));
    
        // Wall Diss. 99
        /*
                      xi2 ^
                          |   
                        1 |   +-----------
                          |  / 
                          | /
                          |/
                          +--------------> Re2
                              1
        */
        const double xi_tau_c = DMIN(re2,1.0);
        tau[2] = vel_norm * hk * 0.5 * xi_tau_c /timefac;
      
    }
    else
    {// stabilization parameters for stationary case
    
        // compute tau_Mu    
        const double re_tau_mu = mk * vel_norm * strle / (2.0 * visc);   /* convective : viscous forces */
        const double xi_tau_mu = DMAX(re_tau_mu, 1.0);
        tau[0] = (DSQR(strle)*mk)/(4.0*visc*xi_tau_mu);
 
        // compute tau_Mp
        const double re_tau_mp = mk * vel_norm * hk / (2.0 * visc);      /* convective : viscous forces */
        const double xi_tau_mp = DMAX(re_tau_mp,1.0);
        tau[1] = (DSQR(hk)*mk)/(4.0*visc*xi_tau_mp);    

        // compute tau_C
        const double xi_tau_c = DMIN(re_tau_mp, 1.0);
        tau[2] = 0.5*vel_norm*hk*xi_tau_c;
    }
    return tau;
}


// check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are necessary
bool DRT::Elements::XFluid3::is_higher_order_element(
              const DRT::Element::DiscretizationType  distype) const
{
    bool hoel = true;
    switch (distype)
    {
    case hex8: case hex20: case hex27: case tet10:
        hoel = true;
        break;
    case tet4:
        hoel = false;
        break;
    default:
        dserror("distype unknown!");
    }
    return hoel;
}


//void DRT::Elements::XFluid3::integrateTerm(
//    EnrEleMatrix&        mat_ele,
//    const double****     shp,
//    const int            numnode,
//    const double         value,
//    const int            adf,
//    const int            shpid_a,
//    const int            bdf,
//    const int            shpid_b
//    )
//{
//    const int numenr = 1;
//    for (int aenr=0;aenr<numenr;aenr++)
//        for (int benr=0;benr<numenr;benr++)
//            for (int anode=0;anode<numnode;anode++)
//                for (int bnode=0;bnode<numnode;bnode++)
//                     mat_ele[aenr][adf][anode][benr][bdf][bnode] 
//                            += shp[aenr][shpid_a][anode][adf] * value * shp[benr][shpid_b][bnode][bdf];
//    return;
//}

//tuple<int, int, int> DRT::Elements::XFluid3::assign_dof_ele(
//    EnrEleMatrix&        mat_ele,
//    const double****     shp,
//    const int            numnode,
//    const double         value,
//    const int            adf,
//    const int            shpid_a,
//    const int            bdf,
//    const int            shpid_b
//    )
//{
//    for (int ienr=0;ienr<numenr;ienr++)
//    {
//        for (int inode=0;inode<numnode;inode++)
//        {
//            for (int idf=0;idf<NDF_;idf++)
//            {
//                if (eleenr[ienr])
//                     mat_ele[aenr][adf][anode][benr][bdf][bnode] 
//                            += shp[aenr][shpid_a][anode][adf] * value * shp[benr][shpid_b][bnode][bdf];
//    return;
//}


// converts a string into an Action for this element
DRT::Elements::XFluid3::ActionType DRT::Elements::XFluid3::convertStringToActionType(
              const string& action) const
{
    
    dsassert(action != "none", "No action supplied");
    
    DRT::Elements::XFluid3::ActionType act = XFluid3::none;
    // get the action required
    
    if (action == "calc_fluid_systemmat_and_residual")      
        act = XFluid3::calc_fluid_systemmat_and_residual;
    else if (action == "calc_fluid_beltrami_error")      
        act = XFluid3::calc_fluid_beltrami_error;
    else if (action == "calc_Shapefunction")
        act = XFluid3::calc_Shapefunction;
    else if (action == "calc_ShapeDeriv1")
        act = XFluid3::calc_ShapeDeriv1;
    else if (action == "calc_ShapeDeriv2")
        act = XFluid3::calc_ShapeDeriv2;
    else 
        dserror("Unknown type of action for XFluid3");
    return act;
}


/*----------------------------------------------------------------------*
 |  init the element (public)                                mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::Elements::XFluid3Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
