/*!----------------------------------------------------------------------
\file fluid3_surface_evaluate.cpp
\brief

Integrate a Surface Neumann boundary condition on a given boundary
element (tri or quad)

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fluid3_xfem.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::XFluid3Surface::Evaluate(	ParameterList& params,
                                    			DRT::Discretization&      discretization,
                                    			vector<int>&              lm,
                                   	 			Epetra_SerialDenseMatrix& elemat1,
                                    			Epetra_SerialDenseMatrix& elemat2,
                                    			Epetra_SerialDenseVector& elevec1,
                                    			Epetra_SerialDenseVector& elevec2,
                                    			Epetra_SerialDenseVector& elevec3)
{
	DRT::Elements::XFluid3Surface::ActionType act = XFluid3Surface::none;
	string action = params.get<string>("action","none");
	if (action == "none") dserror("No action supplied");
	else if (action == "calc_ShapefunctDeriv1Deriv2")
  		act = XFluid3Surface::calc_ShapefunctDeriv1Deriv2;
  	else dserror("Unknown type of action for XFluid3_Surface");
  	
  switch(act)
  {
		case calc_ShapefunctDeriv1Deriv2:
      {
      	// functions, deriv, iel, r, s,               deriv2 is mssing !!!
      	f3_shapefunction_for_surface(elevec1,elemat1,lm[0],elevec2[0],elevec2[1]);
      }
      default:
        dserror("Unknown type of action for XFluid3_Surface");
  } // end of switch(act)

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::XFluid3Surface::EvaluateNeumann(
                                           ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{  
  // there are 3 velocities and 1 pressure 
  const int numdf = 4;

  const double thsl = params.get("time constant for integration",0.0);
  
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  // get the factor for the timecurve
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::TimeCurveManager::Instance().Curve(curvenum).f(time);

  // get values and switches from the condition
  vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  vector<double>* val   = condition.Get<vector<double> >("val"  );

  // set the number of gausspoints
  int nir   = 0;
  int nis   = 0;

  // set number of nodes
  int iel   = this->NumNode();
  
  switch(this->Shape())
  {
      // the surface element is rectangular
      case quad4:
        nir  = 2; // precision 3
        nis  = 2;
        break;
      case quad8:
      case quad9:
        nir  = 3; // precision 5
        nis  = 3;
        break;
      // the surface element is triangular
      case tri3 :
        nir  = 3; /* precision 2               */
        nis  = 1; /* second loop is 'inactive' */
        break;
      case tri6:
        nir  = 6; /* precision 4               */ 
        nis  = 1; /* second loop is 'inactive' */
        break;
      default: dserror("parent element type unknown!\n");
  }


  vector<double>            gaussweight(nir*nis  );  
  Epetra_SerialDenseMatrix  gausscoord (nir*nis,2);


  switch(parent_->Shape())
  {
      case hex8 :
        gaussweight[0]  =  1.0;
        gaussweight[1]  =  1.0;
        gaussweight[2]  =  1.0;
        gaussweight[3]  =  1.0;
        
        gausscoord(0,0) = -0.5773502691896;
        gausscoord(0,1) = -0.5773502691896;
        gausscoord(1,0) =  0.5773502691896;
        gausscoord(1,1) = -0.5773502691896;
        gausscoord(2,0) = -0.5773502691896;
        gausscoord(2,1) =  0.5773502691896;
        gausscoord(3,0) =  0.5773502691896;
        gausscoord(3,1) =  0.5773502691896;
        break;
      case hex20:
      case hex27:
        gaussweight[0]  =  0.5555555555556*0.5555555555556;
        gaussweight[1]  =  0.8888888888889*0.5555555555556;
        gaussweight[2]  =  0.5555555555556*0.5555555555556;
        gaussweight[3]  =  0.5555555555556*0.8888888888889;
        gaussweight[4]  =  0.8888888888889*0.8888888888889;
        gaussweight[5]  =  0.5555555555556*0.8888888888889;
        gaussweight[6]  =  0.5555555555556*0.5555555555556;
        gaussweight[7]  =  0.8888888888889*0.5555555555556;
        gaussweight[8]  =  0.5555555555556*0.5555555555556;
        
        gausscoord(0,0) = -0.7745966692415;
        gausscoord(0,1) = -0.7745966692415;
        gausscoord(1,0) =  0.0;
        gausscoord(1,1) = -0.7745966692415;
        gausscoord(2,0) =  0.7745966692415;
        gausscoord(2,1) = -0.7745966692415;
        gausscoord(3,0) = -0.7745966692415;
        gausscoord(3,1) =  0.0; 
        gausscoord(4,0) =  0.0; 
        gausscoord(4,1) =  0.0;
        gausscoord(5,0) =  0.7745966692415; 
        gausscoord(5,1) =  0.0;
        gausscoord(6,0) = -0.7745966692415; 
        gausscoord(6,1) =  0.7745966692415; 
        gausscoord(7,0) =  0.0;  
        gausscoord(7,1) =  0.7745966692415; 
        gausscoord(8,0) =  0.7745966692415;
        gausscoord(8,1) =  0.7745966692415; 
        break;
      case tet4 :
        gaussweight[0]  = 1.0/6.0 ;
        gaussweight[1]  = 1.0/6.0 ;
        gaussweight[2]  = 1.0/6.0 ;

        gausscoord(0,0) = 0.5;
        gausscoord(0,1) = 0.0;
        gausscoord(1,0) = 0.5;
        gausscoord(1,1) = 0.5;
        gausscoord(2,0) = 0.0;
        gausscoord(2,1) = 0.5;
        break;
      case tet10:
        gaussweight[0]  = 0.0549758718277;
        gaussweight[1]  = 0.0549758718277;
        gaussweight[2]  = 0.0549758718277;
        gaussweight[3]  = 0.1116907948390;
        gaussweight[4]  = 0.1116907948390;
        gaussweight[5]  = 0.1116907948390;

        gausscoord(0,0) = 0.0915762135098;
        gausscoord(0,1) = 0.0915762135098;
        gausscoord(1,0) = 0.8168475729805;
        gausscoord(1,1) = 0.0915762135098;
        gausscoord(2,0) = 0.0915762135098;
        gausscoord(2,1) = 0.8168475729805;
        gausscoord(3,0) = 0.4459484909160;
        gausscoord(3,1) = 0.1081030181681;
        gausscoord(4,0) = 0.4459484909160;
        gausscoord(4,1) = 0.4459484909160;
        gausscoord(5,0) = 0.1081030181681; 
        gausscoord(5,1) = 0.4459484909160; 
        break;
      default: dserror("parent element type unknown!\n");
  }


  // allocate vector for shape functions and matrix for derivatives
  Epetra_SerialDenseVector  funct(iel);
  Epetra_SerialDenseMatrix  deriv(2,iel);

  // node coordinates
  Epetra_SerialDenseMatrix      xyze(3,iel);

  // the metric tensor and the area of an infintesimal surface element
  Epetra_SerialDenseMatrix  metrictensor(2,2);
  double                        drs;

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=this->Nodes()[i]->X()[0];
    xyze(1,i)=this->Nodes()[i]->X()[1];
    xyze(2,i)=this->Nodes()[i]->X()[2];
  }

  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  
  for (int gpid=0;gpid<nir*nis;gpid++)
  {
    double e[2];
    e[0] = gausscoord(gpid,0);
    e[1] = gausscoord(gpid,1);
    
    // get shape functions and derivatives in the plane of the element
    f3_shapefunction_for_surface(funct,deriv,iel,e[0],e[1]);

    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration
    f3_metric_tensor_for_surface(xyze,deriv,metrictensor,&drs);

    // values are multiplied by the product from inf. area element,
    // the gauss weight, the timecurve factor and the constant
    // belonging to the time integration algorithm (theta*dt for
    // one step theta, 2/3 for bdf with dt const.)
    double fac;
    fac  = gaussweight[gpid] * drs * curvefac * thsl;

    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<3;dim++)
      {
        elevec1[node*numdf+dim]+=
          funct[node] * (*onoff)[dim] * (*val)[dim] * fac;
      }
    }

  } /* end of loop over integration points gpid */

  return 0;
}
/*----------------------------------------------------------------------*
get shape function of surface (private) gammi                     04/07 

In this routine the shape functions (always) and their natural first
derivatives with respect to r/s are evaluated for
R E C T A N G L E S or T R I A N G L E S

Numbering of the nodes:


                    ^ s
                    |
              1     |4    0
              o-----o-----o
              |           |
              |           |7
            5 o     o     o -----> r
              |     8     |
              |           |
              o-----o-----o
              2     6     3



                    ^ s
                    |
                    |
                   2|
                    o
                    |\
                    | \ 
                   5o  o4
                    |   \
                    |    \
                    o--o--o -----> r
               0   3   1


                                                            
  
 \param   funct    vector<double>&             (o)    shape functions
 \param   deriv    Epetra_SerialDenseMatrix&   (o)    1st natural deriv.
                                                      of shape funct.
 \param   deriv2   Epetra_SerialDenseMatrix&   (o)    2nd natural deriv.
                                                      of shape funct.
 \param   iel      const int&                  (i)    number of nodes
 \param   r        DOUBLE&                     (i)    coordinate
 \param   s        DOUBLE&                     (i)    coordinate

 *----------------------------------------------------------------------*/



void DRT::Elements::XFluid3Surface::f3_shapefunction_for_surface(	
											Epetra_SerialDenseVector&	funct ,
  											Epetra_SerialDenseMatrix& 	deriv ,
  											const int&                 		iel   ,
  											const double&              		r     ,
  											const double&              		s		)
{
	double Q12=0.50;
  	double Q14=0.25;
  
  	/*------------------------------- selection of polynomial interpolation */
  	switch (iel)
  	{
    	case 4: /* LINEAR shape functions for quad4 and their natural derivatives ----*/
    	{    
    		/*--------------------------------------------- form basic values */
    		double rp=1.0+r;
    		double rm=1.0-r;
    		double sp=1.0+s;
    		double sm=1.0-s;
    
    		funct[0]=Q14*rp*sp;
    		funct[1]=Q14*rm*sp;
   		funct[2]=Q14*rm*sm;
    		funct[3]=Q14*rp*sm;
    
    		deriv(0,0)= Q14*sp;
     		deriv(1,0)= Q14*rp;
        
    		deriv(0,1)=-Q14*sp;
     		deriv(1,1)= Q14*rm;
        
    		deriv(0,2)=-Q14*sm;
     		deriv(1,2)=-Q14*rm;
        
     		deriv(0,3)= Q14*sm;
    		deriv(1,3)=-Q14*rp;
    		break;
    	}
    	case 8: /* QUADRATIC shape functions for quadrilaterals without
           central node and their natural derivatives (serendipity) */
    	{
    		double rp=1.0+r;
    		double rm=1.0-r;
    		double sp=1.0+s;
    		double sm=1.0-s;
    		double r2=1.0-r*r;
    		double s2=1.0-s*s;

			funct[0]=Q14*rp*sp-Q12*(funct[4]+funct[7]);
    		funct[1]=Q14*rm*sp-Q12*(funct[4]+funct[5]);
    		funct[2]=Q14*rm*sm-Q12*(funct[5]+funct[6]);
    		funct[3]=Q14*rp*sm-Q12*(funct[6]+funct[7]);
    		funct[4]=Q12*r2*sp;
    		funct[5]=Q12*rm*s2;
    		funct[6]=Q12*r2*sm;
    		funct[7]=Q12*rp*s2;
    	
      	deriv(0,0)= Q14*sp;
      	deriv(1,0)= Q14*rp;
        
      	deriv(0,1)=-Q14*sp;
      	deriv(1,1)= Q14*rm;
        
      	deriv(0,2)=-Q14*sm;
      	deriv(1,2)=-Q14*rm;
        
      	deriv(0,3)= Q14*sm;
      	deriv(1,3)=-Q14*rp;
        
      	deriv(0,4)=-1.0*r*sp;
      	deriv(1,4)= Q12*r2;
        
      	deriv(0,5)=-Q12*  s2;
      	deriv(1,5)=-1.0*rm*s;
        
      	deriv(0,6)=-1.0*r*sm;
      	deriv(1,6)=-Q12*r2;
        
      	deriv(0,7)= Q12*s2;
      	deriv(1,7)=-1.0*rp*s;
        
      	deriv(0,0)-= Q12*(deriv(0,4)+deriv(0,7));
      	deriv(1,0)-= Q12*(deriv(1,4)+deriv(1,7));

    		for(int i=1;i<4;i++)
   		{
      		int ii=i+3;
         	deriv(0,i) -= Q12*(deriv(0,ii)+deriv(0,ii+1));
         	deriv(1,i) -= Q12*(deriv(1,ii)+deriv(1,ii+1));
			} /* end loop over i */
    		break;
    	}
    	case 9: /* full QUADRATIC shape functions for quadrilaterals with
                         central node and their natural derivatives */
    	{
			/*--------------------------------------------------- form basic values */
			double rp=1.0+r;
	   	double rm=1.0-r;
	   	double sp=1.0+s;
	   	double sm=1.0-s;
	   	double r2=1.0-r*r;
	   	double s2=1.0-s*s;
      	double rh=Q12*r;
      	double sh=Q12*s;
      	double rs=rh*sh;
	   	double rhp=r+Q12;
	   	double rhm=r-Q12;
	   	double shp=s+Q12;
	   	double shm=s-Q12;
	
	   	funct[0]= rs*rp*sp;
	   	funct[1]=-rs*rm*sp;
	   	funct[2]= rs*rm*sm;
	   	funct[3]=-rs*rp*sm;
		   funct[4]= sh*sp*r2;
		   funct[5]=-rh*rm*s2;
		   funct[6]=-sh*sm*r2;
		   funct[7]= rh*rp*s2;
		   funct[8]= r2*s2;
	
			deriv(0,0)= rhp*sh*sp;
		   deriv(1,0)= shp*rh*rp;
		     
		   deriv(0,1)= rhm*sh*sp;
		   deriv(1,1)=-shp*rh*rm;
		     
		   deriv(0,2)=-rhm*sh*sm;
		   deriv(1,2)=-shm*rh*rm;
		     
		   deriv(0,3)=-rhp*sh*sm;
		   deriv(1,3)= shm*rh*rp;
		     
		   deriv(0,4)=-2.0*r*sh*sp;
		   deriv(1,4)= shp*r2;
		     
	      deriv(0,5)= rhm*s2;
	      deriv(1,5)= 2.0*s*rh*rm;
	     
	      deriv(0,6)= 2.0*r*sh*sm;
	      deriv(1,6)= shm*r2;
	     
	      deriv(0,7)= rhp*s2;
	      deriv(1,7)=-2.0*s*rh*rp;      
			
	      deriv(0,8)=-2.0*r*s2;
	      deriv(1,8)=-2.0*s*r2;
	    break;
    	}
    	case 3: /* LINEAR shape functions for triangles and their natural derivatives -----*/
    	{
			/*------------------------------------------- form basic values */
    		funct[0]=1.0-r-s;
    		funct[1]=r;
    		funct[2]=s;
    
      	deriv(0,0)=-1.0;
      	deriv(1,0)=-1.0;
      	deriv(0,1)= 1.0;
      	deriv(1,1)= 0.0;
      	deriv(0,2)= 0.0;
      	deriv(1,2)= 1.0;
    		break;
    	}
    	case 6: /* QUADRATIC shape functions for triangles and their natural derivatives -*/
    	{
      	  /*------------------------------------------- form basic values */
	    	double rr=r*r;
	    	double ss=s*s;
	    	double rs=r*s;
	
	    	funct[0]=(1.0-2.0*r-2.0*s)*(1.0-r-s);
	    	funct[1]=2.0*rr-r;
		    funct[2]=2.0*ss-s;
		    funct[3]=4.0*(r-rr-rs);
		    funct[4]=4.0*rs;
		    funct[5]=4.0*(s-rs-ss);
	    
	       deriv(0,0)=-3.0+4.0*(r+s);
	       deriv(1,0)= deriv(0,0);
	        
	       deriv(0,1)= 4.0*r-1.0;
	       deriv(1,1)= 0.0;
	        
	       deriv(0,2)= 0.0;
	       deriv(1,2)= 4.0*s-1.0;
	        
	       deriv(0,3)= 4.0*(1.0-2.0*r-s);
	       deriv(1,3)=-4.0*r;
	        
	       deriv(0,4)= 4.0*s;
	       deriv(1,4)= 4.0*r;
	        
	       deriv(0,5)=-4.0*s;
	       deriv(1,5)= 4.0*(1.0-r-2.0*s);
	    	break;
	    }
    	/*------------------------------------------------------------------*/
    	default:
    		dserror("distyp unknown\n");
 	} /* end switch(iel) */
 
  	return;
}

/* compute kovariant metric tensor G for fluid element        gammi 04/07

                        +-       -+
                        | g11 g12 |
                    G = |         |
                        | g12 g22 |
                        +-       -+

 where (o denotes the inner product, xyz a vector)


                            dxyz   dxyz
                    g11 =   ---- o ----
                             dr     dr
                    
                            dxyz   dxyz
                    g12 =   ---- o ----
                             dr     ds
                    
                            dxyz   dxyz
                    g22 =   ---- o ----
                             ds     ds
                    
                    
 and the square root of the first fundamental form


                          +--------------+
                         /               |
           sqrtdetg =   /  g11*g22-g12^2
                      \/

 they are needed for the integration over the surface element
               
*/
void  DRT::Elements::XFluid3Surface::f3_metric_tensor_for_surface(
  const Epetra_SerialDenseMatrix  xyze,
  const Epetra_SerialDenseMatrix  deriv,
  Epetra_SerialDenseMatrix&       metrictensor,
  double                         *sqrtdetg)
{
  /* 
  |                                              0 1 2
  |                                             +-+-+-+ 
  |       0 1 2              0...iel-1          | | | | 0
  |      +-+-+-+             +-+-+-+-+          +-+-+-+ 
  |      | | | | 1           | | | | | 0        | | | | .
  |      +-+-+-+       =     +-+-+-+-+       *  +-+-+-+ .
  |      | | | | 2           | | | | | 1        | | | | .
  |      +-+-+-+             +-+-+-+-+          +-+-+-+ 
  |                                             | | | | iel-1
  |                                     +-+-+-+ 
  |                                      
  |       dxyzdrs             deriv              xyze^T
  |                                        
  |                                        
  |                                     +-            -+
  |                                 | dx   dy   dz |
  |                                 | --   --   -- |
  |                         | dr   dr   dr |
  |     yields               dxyzdrs =  |              |
  |                                 | dx   dy   dz |
  |                                 | --   --   -- |
  |                         | ds   ds   ds |
  |                                     +-            -+
  |
  */
  Epetra_SerialDenseMatrix dxyzdrs (2,3);
  
  dxyzdrs.Multiply('N','T',1.0,deriv,xyze,0.0);

  /*
  |                                        
  |      +-           -+    +-            -+   +-            -+ T
  |      |             |    | dx   dy   dz |   | dx   dy   dz |
  |      |  g11   g12  |    | --   --   -- |   | --   --   -- |
  |      |             |    | dr   dr   dr |   | dr   dr   dr |
  |      |             |  = |              | * |              | 
  |      |             |    | dx   dy   dz |   | dx   dy   dz |
  |      |  g21   g22  |    | --   --   -- |   | --   --   -- |
  |      |             |    | ds   ds   ds |   | ds   ds   ds |
  |      +-           -+    +-            -+   +-            -+
  |
  | the calculation of g21 is redundant since g21=g12
  */
  metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,0.0);

/*
                          +--------------+
                         /               |
           sqrtdetg =   /  g11*g22-g12^2
                      \/
*/

  sqrtdetg[0]= sqrt(metrictensor(0,0)*metrictensor(1,1)
                    -
                    metrictensor(0,1)*metrictensor(1,0));

  return;
}

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
