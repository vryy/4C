/*!----------------------------------------------------------------------
\file fluid2_line_evaluate.cpp
\brief

<pre>
Maintainer: Peter Gmanitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fluid2.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"

using namespace DRT::Utils;


/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)     gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::Fluid2Line::EvaluateNeumann(
    ParameterList& params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1)
{  
  
  // there are 2 velocities and 1 pressure
  const int numdf = 3;

  const double thsl = params.get("time constant for integration",0.0);

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
  // (assumed to be constant on element boundary)
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );
  

    // set number of nodes
  const int iel   = this->NumNode();
    
  const DiscretizationType distype = this->Shape();
  
  // gaussian points
  const GaussRule1D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints1D  intpoints = getIntegrationPoints1D(gaussrule);
  
 
  // allocate vector for shape functions and for derivatives
  Epetra_SerialDenseVector   funct(iel);
  Epetra_SerialDenseMatrix    deriv(iel,1);

  // node coordinates
  Epetra_SerialDenseMatrix xye(2,iel);

  // the length of an infintesimal line element
  double dr=0;

  // get node coordinates
  for(int i=0;i<iel;++i)
  {
    xye(0,i)=this->Nodes()[i]->X()[0];
    xye(1,i)=this->Nodes()[i]->X()[1];
  }
  
  
  // loop over integration points
  for (int gpid=0;gpid<intpoints.nquad;gpid++)
  {
    const double e1 = intpoints.qxg[gpid];
    // get shape functions and derivatives in the line
    shape_function_1D(funct,e1,distype);
    shape_function_1D_deriv1(deriv,e1,distype);

    // compute infintesimal line element dr for integration along the line
    f2_substitution(xye,deriv,dr,iel);
  
    // values are multiplied by the product from inf. area element,
    // the gauss weight, the timecurve factor and the constant
    // belonging to the time integration algorithm (theta*dt for
    // one step theta, 2/3 for bdf with dt const.)
    
    double fac;
    fac =intpoints.qwgt[gpid] *dr* curvefac * thsl;
     for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<3;dim++)
      {
        elevec1[node*numdf+dim]+=
          funct[node] * (*onoff)[dim] * (*val)[dim] * fac;
      }
    }
  }//end of loop over integrationen points

  //dserror("Line Neumann condition not yet implemented for Fluid2");
  return 0;
}

GaussRule1D DRT::Elements::Fluid2Line::getOptimalGaussrule(const DiscretizationType& distype)
{
  GaussRule1D rule;
  switch (distype)
    {
    case line2:
      rule = intrule_line_2point;
      break;
    case line3:
      rule = intrule_line_3point;
      break;
    default: 
    dserror("unknown number of nodes for gaussrule initialization");
    }
  return rule;
}


void  DRT::Elements::Fluid2Line::f2_substitution(
  const Epetra_SerialDenseMatrix  xye,
  const Epetra_SerialDenseMatrix  deriv,
  double&               dr,
  const int iel)
{
  // compute derivative of parametrization
  Epetra_SerialDenseVector der_par (iel);
  der_par.Multiply('N','N',1.0,xye,deriv,0.0);
  dr=der_par.Norm2();
  return;

}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID2
