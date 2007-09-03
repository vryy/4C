/*!----------------------------------------------------------------------
\file so_weg6_surface.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOH8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "so_weg6.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"

extern "C"
{
#include "../headers/standardtypes.h"
}
#include "../drt_lib/dstrc.H"

using namespace DRT::Utils;

/*----------------------------------------------------------------------*
 * Integrate a Surface Neumann boundary condition (public)     maf 04/07*
 * ---------------------------------------------------------------------*/
int DRT::Elements::Sow6Surface::EvaluateNeumann(ParameterList&           params,
                                                DRT::Discretization&     discretization,
                                                DRT::Condition&          condition,
                                                vector<int>&             lm,
                                                Epetra_SerialDenseVector& elevec1)
{
  // get type of condition
  enum LoadType
  {
    neum_none,
    neum_live,
    neum_orthopressure,
    neum_consthydro_z,
    neum_increhydro_z,
    neum_live_FSI,
    neum_opres_FSI
  };
  LoadType ltype;
  const string* type = condition.Get<string>("type");
  if      (*type == "neum_live")          ltype = neum_live;
  //else if (*type == "neum_live_FSI")      ltype = neum_live_FSI;
  else dserror("Unknown type of SurfaceNeumann condition");

  // get values and switches from the condition
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );

  /*
  **    TIME CURVE BUSINESS
  */
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
  // **

  // element geometry update
  const DiscretizationType distype = this->Shape();
  const int numnod = this->NumNode();
  
  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
  
  Epetra_SerialDenseMatrix xsrefe(numnod,NUMDIM_WEG6);  // material coord. of element
  Epetra_SerialDenseMatrix xscurr(numnod,NUMDIM_WEG6);  // material coord. of element
  for (int i=0; i<numnod; ++i){
    xsrefe(i,0) = Nodes()[i]->X()[0];
    xsrefe(i,1) = Nodes()[i]->X()[1];
    xsrefe(i,2) = Nodes()[i]->X()[2];
    
    xscurr(i,0) = xsrefe(i,0) + mydisp[i*NODDOF_WEG6+0];
    xscurr(i,1) = xsrefe(i,1) + mydisp[i*NODDOF_WEG6+1];
    xscurr(i,2) = xsrefe(i,2) + mydisp[i*NODDOF_WEG6+2];
  }

  switch(distype){
  case quad4:
    gaussrule = intrule_quad_4point;
    break;
  case tri3:
    gaussrule = intrule_tri_3point;
    break;
  default: 
      dserror("shape type unknown!\n");
  }
  
  // allocate vector for shape functions and matrix for derivatives
  Epetra_SerialDenseVector  funct       (numnode);
  Epetra_SerialDenseMatrix  deriv       (2,numnode);

  // the metric tensor and the area of an infintesimal surface element
  Epetra_SerialDenseMatrix  metrictensor(2,2);
  double                        drs;

  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  const IntegrationPoints2D  intpoints = getIntegrationPoints2D(gaussrule);
  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    shape_function_2D(funct, e0, e1, distype);
    shape_function_2D_deriv1(deriv, e0, e1, distype);

    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration
    sow6_surface_integ(&drs,xsrefe,deriv);

    // values are multiplied by the product from inf. area element,
    // the gauss weight, the timecurve factor
    const double fac = intpoints.qwgt[gpid] * drs * curvefac;

    for (int node=0; node < numnod; ++node)
    {
      for(int dim=0 ; dim<NUMDIM_WEG6; dim++)
      {
        elevec1[node*numdf+dim]+=
          funct[node] * (*onoff)[dim] * (*val)[dim] * fac;
      }
    }

  } /* end of loop over integration points gpid */
  
  return 0;
}

/*----------------------------------------------------------------------*
 * Evaluate sqrt of determinant of metric at gp (private)      maf 05/07*
 * ---------------------------------------------------------------------*/
void DRT::Elements::Sow6Surface::sow6_surface_integ(
      double* sqrtdetg,                      // (o) pointer to sqrt of det(g)
      const Epetra_SerialDenseMatrix* xs,    // (i) element coords
      const Epetra_SerialDenseMatrix deriv)  // (i) shape funct derivs
{

  // compute dXYZ / drs
  Epetra_SerialDenseMatrix dxyzdrs(2,3);
  dxyzdrs.Multiply('T','N',1.0,deriv,(*xs),1.0);

  /* compute covariant metric tensor G for surface element
  **                        | g11   g12 |
  **                    G = |           |
  **                        | g12   g22 |
  ** where (o denotes the inner product, xyz a vector)
  **
  **       dXYZ   dXYZ          dXYZ   dXYZ          dXYZ   dXYZ
  ** g11 = ---- o ----    g12 = ---- o ----    g22 = ---- o ----
  **        dr     dr            dr     ds            ds     ds
  */
  Epetra_SerialDenseMatrix metrictensor(2,2);
  metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,1.0);
  (*sqrtdetg) = sqrt( metrictensor(0,0)*metrictensor(1,1)
                     -metrictensor(0,1)*metrictensor(1,0));
                       
  return;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOH8
