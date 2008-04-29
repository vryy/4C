/*!----------------------------------------------------------------------
\file so_surface_evaluate.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_surface.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"

/*----------------------------------------------------------------------*
 * Integrate a Surface Neumann boundary condition (public)     gee 04/08|
 * ---------------------------------------------------------------------*/
int DRT::ELEMENTS::StructuralSurface::EvaluateNeumann(ParameterList&           params,
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
    neum_orthopressure
  };
  // spatial or material configuration depends on the type of load
  enum Configuration
  {
    config_none,
    config_material,
    config_spatial,
    config_both
  };
  
  Configuration config = config_none;
  LoadType ltype       = neum_none;
  const string* type = condition.Get<string>("type");
  if      (*type == "neum_live")          
  {
    ltype = neum_live;
    config = config_material;
  }
  else if (*type == "neum_orthopressure") 
  {
    ltype = neum_orthopressure;
    config = config_spatial;
  }
  else 
  {
    dserror("Unknown type of SurfaceNeumann condition");
  }

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
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

  // element geometry update
  const int numnode = NumNode();
  const int numdf=3;
  LINALG::SerialDenseMatrix x(numnode,3);
  LINALG::SerialDenseMatrix xc;
  switch (config)
  {
  case config_material:
    MaterialConfiguration(x);
  break;
  case config_spatial:
  {
    RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    if (disp==null) dserror("Cannot get state vector 'displacement'");
    vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    SpatialConfiguration(x,mydisp);
  }
  break;
  case config_both:
  {
    xc.LightShape(numnode,3);
    RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    if (disp==null) dserror("Cannot get state vector 'displacement'");
    vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    MaterialConfiguration(x);
    SpatialConfiguration(xc,mydisp);
  }
  break;
  default: dserror("Unknown case of frame");
  break;
  }
  
  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseVector  funct(numnode);
  LINALG::SerialDenseMatrix  deriv(2,numnode);

  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  const DRT::UTILS::IntegrationPoints2D  intpoints = 
                                       getIntegrationPoints2D(gaussrule_);
  for (int gp=0; gp<intpoints.nquad; gp++)
  {
    const double e0 = intpoints.qxg[gp][0];
    const double e1 = intpoints.qxg[gp][1];

    // get shape functions and derivatives in the plane of the element
    DRT::UTILS::shape_function_2D(funct,e0,e1,Shape());
    DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

    switch(ltype)
    {
    case neum_live:
    {
      LINALG::SerialDenseMatrix dxyzdrs(2,3);
      dxyzdrs.Multiply('N','N',1.0,deriv,x,0.0);
      LINALG::SerialDenseMatrix  metrictensor(2,2);
      metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,0.0);
      const double detA = sqrt( metrictensor(0,0)*metrictensor(1,1)
                                -metrictensor(0,1)*metrictensor(1,0));
      const double fac = intpoints.qwgt[gp] * detA * curvefac;
      for (int node=0; node < numnode; ++node)
        for(int dim=0 ; dim<3; ++dim)
          elevec1[node*numdf+dim]+= funct[node] * (*onoff)[dim] * (*val)[dim] * fac;
    }
    break;
    case neum_orthopressure:
    {
     if ((*onoff)[0] != 1) dserror("orthopressure on 1st dof only!");
      for (int checkdof = 1; checkdof < 3; ++checkdof) {
        if ((*onoff)[checkdof] != 0) dserror("orthopressure on 1st dof only!");
      }
      double ortho_value = (*val)[0];
      if (!ortho_value) dserror("no orthopressure value given!");
      vector<double> normal(3);
      double detA;
      SurfaceIntegration(detA,normal,x,deriv);
      const double fac = intpoints.qwgt[gp] * curvefac * ortho_value * (-1.0);
      for (int node=0; node < numnode; ++node)
        for(int dim=0 ; dim<3; dim++)
          elevec1[node*numdf+dim]+= funct[node] * normal[dim] * fac;
    }
    break;
    default:
      dserror("Unknown type of SurfaceNeumann load");
    break;
    }

  } /* end of loop over integration points gp */
  
  return 0;
}

/*----------------------------------------------------------------------*
 * Evaluate sqrt of determinant of metric at gp (private)      gee 04/08|
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::SurfaceIntegration(double& detA,                      
                                                          vector<double>& normal,                  
                                                          const Epetra_SerialDenseMatrix& x,    
                                                          const Epetra_SerialDenseMatrix& deriv)  
{

  // compute dXYZ / drs
  LINALG::SerialDenseMatrix dxyzdrs(2,3);
  dxyzdrs.Multiply('N','N',1.0,deriv,x,0.0);

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
  LINALG::SerialDenseMatrix metrictensor(2,2);
  metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,0.0);
  detA = sqrt( metrictensor(0,0)*metrictensor(1,1)-metrictensor(0,1)*metrictensor(1,0) );
  normal[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
  normal[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
  normal[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);
                       
  return;
}


#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOLID3
