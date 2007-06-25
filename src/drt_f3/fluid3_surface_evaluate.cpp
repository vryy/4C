/*!----------------------------------------------------------------------
\file fluid3_surface_evaluate.cpp
\brief

Integrate a Surface Neumann boundary condition on a given boundary
element (tri or quad)

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fluid3.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"



/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::Fluid3Surface::EvaluateNeumann(
                                           ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  // there are 3 velocities and 1 pressure
  const int numdf = 4;

  const double thsl = params.get("time constant for integration",0.0);
  
  const DiscretizationType distype = this->Shape();

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(time);

  // get values and switches from the condition
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );

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
  Epetra_SerialDenseVector  funct       (iel);
  Epetra_SerialDenseMatrix 	deriv       (2,iel);

  // node coordinates
  Epetra_SerialDenseMatrix      xyze        (3,iel);

  // the metric tensor and the area of an infintesimal surface element
  Epetra_SerialDenseMatrix 	metrictensor(2,2);
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
    DRT::Utils::shape_function_2D(funct,e[0],e[1],distype);
    DRT::Utils::shape_function_2D_deriv1(deriv,e[0],e[1],distype);

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
void  DRT::Elements::Fluid3Surface::f3_metric_tensor_for_surface(
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
  |		     	      	     	        +-+-+-+
  |
  |       dxyzdrs             deriv              xyze^T
  |
  |
  |                                     +-            -+
  |  	   	    	    	        | dx   dy   dz |
  |  	   	    	    	        | --   --   -- |
  | 	   	   	   	        | dr   dr   dr |
  | 	yields               dxyzdrs =  |              |
  |  	   	    	    	        | dx   dy   dz |
  |  	   	    	    	        | --   --   -- |
  | 	   	   	   	        | ds   ds   ds |
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
