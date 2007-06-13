/*!----------------------------------------------------------------------
\file so_hex8_surface.cpp
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

#include "so_hex8.H"
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

/*----------------------------------------------------------------------*
 * Integrate a Surface Neumann boundary condition (public)     maf 04/07*
 * ---------------------------------------------------------------------*/
int DRT::Elements::Soh8Surface::EvaluateNeumann(ParameterList&           params,
                                                DRT::Discretization&     discretization,
                                                DRT::Condition&          condition,
                                                vector<int>&             lm,
                                                Epetra_SerialDenseVector& elevec1)
{
  DSTraceHelper dst("Soh8Surface::EvaluateNeumann");

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
    curvefac = DRT::TimeCurveManager::Instance().Curve(curvenum).f(time);
  // **

  // element geometry
  const int numnod = 4;
  Epetra_SerialDenseMatrix xsrefe(numnod,NUMDIM_SOH8);  // material coord. of element
  for (int i=0; i<numnod; ++i){
    xsrefe(i,0) = Nodes()[i]->X()[0];
    xsrefe(i,1) = Nodes()[i]->X()[1];
    xsrefe(i,2) = Nodes()[i]->X()[2];
  }

  /*
  ** Here, we integrate a 4-node surface with 2x2 Gauss Points
  */
  const int ngp = 4;

  // gauss parameters
  const double gpweight = 1.0;
  const double gploc    = 1.0/sqrt(3.0);
  Epetra_SerialDenseMatrix gpcoord (ngp,2);
  gpcoord(0,0) = - gploc;
  gpcoord(0,1) = - gploc;
  gpcoord(1,0) =   gploc;
  gpcoord(1,1) = - gploc;
  gpcoord(2,0) = - gploc;
  gpcoord(2,1) =   gploc;
  gpcoord(3,0) =   gploc;
  gpcoord(3,1) =   gploc;

  for (int gpid = 0; gpid < 4; ++gpid) {    // loop over intergration points
    // get shape functions and derivatives of element surface
    vector<double> funct(ngp);                // 4 shape function values
    double drs;                               // surface area factor
    soh8_surface_integ(&funct,&drs,&xsrefe,gpcoord(gpid,0),gpcoord(gpid,1));
    double fac = gpweight * drs * curvefac;   // integration factor

    // distribute over element load vector
    for (int nodid=0; nodid < 4; ++nodid) {
      for(int dim=0; dim < NUMDIM_SOH8; ++dim) {
//        int on_off=(*onoff)[dim];
//        double value=(*val)[dim];
//        cout << "id: " << nodid*NUMDOF_SOH8 + dim << "onoff: " << on_off << " value: " << value << endl;
        elevec1[nodid*NUMDIM_SOH8 + dim] += funct[nodid] * (*onoff)[dim] * (*val)[dim] * fac;
      }
    }
  }
//    cout << elevec1 << endl;
  return 0;
}

/*----------------------------------------------------------------------*
 * Evaluate sqrt of determinant of metric at gp (private)      maf 05/07*
 * ---------------------------------------------------------------------*/
void DRT::Elements::Soh8Surface::soh8_surface_integ(
      vector<double>* funct,                 // (o) shape functions
      double* sqrtdetg,                      // (o) pointer to sqrt of det(g)
      const Epetra_SerialDenseMatrix* xsrefe,// (i) material element coords
      const double r,                        // (i) coord in r-direction
      const double s)                        // (i) coord in s-direction
{
  DSTraceHelper dst("Soh8Surface::soh8_surface_metric");

  // shape functions for 4 nodes
  (*funct)[0] = 0.25 * (1.0-r) * (1.0-s);
  (*funct)[1] = 0.25 * (1.0+r) * (1.0-s);
  (*funct)[2] = 0.25 * (1.0+r) * (1.0+s);
  (*funct)[3] = 0.25 * (1.0-r) * (1.0+s);
  // derivatives of 4 shape functions wrt 2 directions
  Epetra_SerialDenseMatrix deriv(4,2);
  deriv(0,0) = -0.25 * (1.0-s);
  deriv(0,1) = -0.25 * (1.0-r);
  deriv(1,0) =  0.25 * (1.0-s);
  deriv(1,1) = -0.25 * (1.0+r);
  deriv(2,0) =  0.25 * (1.0+s);
  deriv(2,1) =  0.25 * (1.0+r);
  deriv(3,0) = -0.25 * (1.0+s);
  deriv(3,1) =  0.25 * (1.0-r);

  // compute dXYZ / drs
  Epetra_SerialDenseMatrix dxyzdrs(2,3);
  dxyzdrs.Multiply('T','N',1.0,deriv,(*xsrefe),1.0);

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
