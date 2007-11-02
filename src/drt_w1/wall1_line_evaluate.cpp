/*!----------------------------------------------------------------------
\file wall1_line_evaluate.cpp
\brief

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#ifdef CCADISCRET

#include "wall1.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"

extern "C"
{
#include "../headers/standardtypes.h"
#include "../wall1/wall1.h"
}
#include "../drt_lib/dstrc.H"


/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)     mgit 03/07|
 *----------------------------------------------------------------------*/

int DRT::Elements::Wall1Line::EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1)
{
  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::Utils::ExtractMyValues(*disp,mydisp,lm);

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

   // set number of nodes
  const int iel   = this->NumNode();
  
  const DiscretizationType distype = this->Shape();
  
  // gaussian points 
  const DRT::Utils::GaussRule1D gaussrule = getOptimalGaussrule(distype); 
  const DRT::Utils::IntegrationPoints1D  intpoints = getIntegrationPoints1D(gaussrule);
  
  // get values and switches from the condition
  const vector<int>*    onoff = condition.Get<vector<int> >("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val");
   
  // allocate vector for shape functions and for derivatives
  Epetra_SerialDenseVector    funct(iel);
  Epetra_SerialDenseMatrix    deriv(iel,1);
  
  // node coordinates
  Epetra_SerialDenseMatrix xye(2,iel);
    
  // get node coordinates
  for(int i=0;i<iel;++i)
  {
    xye(0,i)=this->Nodes()[i]->X()[0]; 
    xye(1,i)=this->Nodes()[i]->X()[1];
  }
   
  // loop over integration points //new
  for (int gpid=0;gpid<intpoints.nquad;gpid++)
  {  
     const double e1 = intpoints.qxg[gpid];
     
  // get shape functions and derivatives in the line
     DRT::Utils::shape_function_1D(funct,e1,distype);
     DRT::Utils::shape_function_1D_deriv1(deriv,e1,distype);
  
  // compute infinitesimal line element dr for integration along the line
     const double dr = w1_substitution(xye,deriv,iel);
    
  // load vector ar
     double ar[2];
  // loop the dofs of a node
  // ar[i] = ar[i] * facr * ds * onoff[i] * val[i]*curvefac
     for (int i=0; i<2; ++i)
     {
        ar[i] = intpoints.qwgt[gpid] * dr * (*onoff)[i]*(*val)[i] * curvefac;
     }

   // add load components
     for (int node=0; node<iel; ++node)
     {	 
       for (int j=0; j<2; ++j)
       {
          elevec1[node*2+j] += funct[node]*ar[j];
       }
     }
  } 
  return 0;
}


DRT::Utils::GaussRule1D DRT::Elements::Wall1Line::getOptimalGaussrule(const DiscretizationType& distype)
{
  DRT::Utils::GaussRule1D rule;
  switch (distype)
    {
    case line2:
      rule = DRT::Utils::intrule_line_2point;
      break;
    case line3:
      rule = DRT::Utils::intrule_line_3point;
      break;
    default: 
    dserror("unknown number of nodes for gaussrule initialization");
    }
  return rule;
}

// determinant of jacobian matrix

double  DRT::Elements::Wall1Line::w1_substitution(const Epetra_SerialDenseMatrix& xye,
                                                  const Epetra_SerialDenseMatrix& deriv,
                                                  const int iel)
{
	  /*
	  |                                            0 1 
	  |                                           +-+-+
	  |       0 1              0...iel-1          | | | 0
	  |      +-+-+             +-+-+-+-+          +-+-+
	  |      | | | 1     =     | | | | | 0        | | | .
	  |      +-+-+             +-+-+-+-+       *  +-+-+ .
	  |                                           | | | .
	  |                                           +-+-+
	  |                                           | | | iel-1
	  |		           	     	                  +-+-+
	  |
	  |       dxyzdrs             deriv^T          xye^T
	  |
	  |
	  |                         +-        -+
	  |  	   	    	    	| dx   dy  |
	  |  	  yields   dxydr =	| --   --  |
	  | 	                    | dr   dr  |
	  |                         +-        -+
	  |
	  */
// compute derivative of parametrization
double dr = 0.0;
Epetra_SerialDenseMatrix der_par (1,2);
int err = der_par.Multiply('T','T',1.0,deriv,xye,0.0);
if (err!=0)
	dserror("Multiply failed");
dr=sqrt(der_par(0,0)*der_par(0,0)+der_par(0,1)*der_par(0,1));
return dr;
}


#endif  // #ifdef CCADISCRET
#endif // #ifdef D_WALL1
