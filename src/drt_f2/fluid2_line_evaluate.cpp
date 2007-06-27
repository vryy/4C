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
  

  // set the number of gausspoints
  int nir   = 0;

  // set number of nodes
  int iel   = this->NumNode();
    
  switch(this->Shape())
  {
    case line2:
      nir = 2;
      break;
    case line3:
      nir = 3;
      break;
    default: 
      dserror("line element type unknown!\n");
  }

  vector<double>            gaussweight(nir);
  vector<double>            gausscoord(nir);
  switch(this->Shape())
  {
    case line2:
      gaussweight[0]=1.0;
      gaussweight[1]=1.0;
      gausscoord[0]=-0.57735026919;
      gausscoord[1]=0.57735026919;
      break;
    case line3:
      gaussweight[0]=0.555555555556;
      gaussweight[2]=0.888888888889;
      gaussweight[1]=0.555555555556;
      gausscoord[0]=-0.774596669241;
      gausscoord[2]=0;
      gausscoord[1]=0.774596669241;
      break;
  default:
     dserror("line element type unknown!\n");
  }
 
  // allocate vector for shape functions and for derivatives
  vector<double>   funct(iel);
  vector<double>   deriv2(iel);  
  Epetra_SerialDenseVector    deriv(iel);

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
  for (int gpid=0;gpid<nir;gpid++)
  {
    // get shape functions and derivatives in the line
    f2_shapefunction_for_line(funct,deriv,iel,gausscoord[gpid]);

    // compute infintesimal line element dr for integration along the line
    f2_substitution(xye,deriv,dr);
  
    // values are multiplied by the product from inf. area element,
    // the gauss weight, the timecurve factor and the constant
    // belonging to the time integration algorithm (theta*dt for
    // one step theta, 2/3 for bdf with dt const.)
    
    double fac;
    fac = gaussweight[gpid]*dr* curvefac * thsl;
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


void DRT::Elements::Fluid2Line::f2_shapefunction_for_line(
  vector<double>&           funct ,
  Epetra_SerialDenseVector&         deriv ,
  const int                 iel   ,
  const double              r
  )
{
  switch (iel)
  {
    case 2:
    {
      funct[0] = (ONE-r)*0.5;
      funct[1] = (ONE+r)*0.5;
      
      deriv[0] = -0.5;
      deriv[1] = 0.5;
      break;
    }
    case 3:
    {
      funct[0] = 0.5*r*(r-ONE);
      funct[1] = 0.5*r*(r+ONE);
      funct[2] = (ONE-r)*(ONE+r);

      deriv[0] = r-0.5;
      deriv[1] = r+0.5;
      deriv[2] = -TWO*r;
      break;
    }
    default:
	dserror("distyp unknown\n");
  }
  return;
}


void  DRT::Elements::Fluid2Line::f2_substitution(
  const Epetra_SerialDenseMatrix  xye,
  const Epetra_SerialDenseVector  deriv,
  double&               dr)
{
  // compute derivative of parametrization
  Epetra_SerialDenseVector der_par (deriv.Length());
  der_par.Multiply('N','N',1.0,xye,deriv,0.0);
  dr=der_par.Norm2();
  return;

}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID2
