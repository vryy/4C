/*!----------------------------------------------------------------------
\file so_hex8_line.cpp
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

#include "so_hex8.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"


/*-----------------------------------------------------------------------*
 * Integrate a Line Neumann boundary condition (public)         maf 04/07*
 * ----------------------------------------------------------------------*/
int DRT::Elements::Soh8Line::EvaluateNeumann(ParameterList&         params,
                                             DRT::Discretization&   discretization,
                                             DRT::Condition&        condition,
                                             vector<int>&           lm,
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
  else dserror("Unknown type of LineNeumann condition");

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
  // **

  // element geometry update
  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
  const int numnod = 2;
  Epetra_SerialDenseMatrix xsrefe(numnod,NUMDIM_SOH8);  // material coord. of element
  Epetra_SerialDenseMatrix xscurr(numnod,NUMDIM_SOH8);  // material coord. of element
  for (int i=0; i<numnod; ++i){
    xsrefe(i,0) = Nodes()[i]->X()[0];
    xsrefe(i,1) = Nodes()[i]->X()[1];
    xsrefe(i,2) = Nodes()[i]->X()[2];

    xscurr(i,0) = xsrefe(i,0) + mydisp[i*NODDOF_SOH8+0];
    xscurr(i,1) = xsrefe(i,1) + mydisp[i*NODDOF_SOH8+1];
    xscurr(i,2) = xsrefe(i,2) + mydisp[i*NODDOF_SOH8+2];
  }

  /*
  ** Here, we integrate a 2-node line with 2 Gauss Points
  */
  const int ngp = 2;

  // gauss parameters
  const double gpweight = 1.0;
  const double gploc    = 1.0/sqrt(3.0);
  double gpcoord[2];
  gpcoord[0] = - gploc;
  gpcoord[1] =   gploc;

  for (int gpid = 0; gpid < 2; ++gpid) {    // loop over intergration points
    // get shape functions and derivatives of element surface
    vector<double> funct(ngp);                // 2 shape function values
    double drs;                               // line length factor
    switch(ltype){
      case neum_live:{            // uniform load on reference configuration
        soh8_line_integ(&funct,&drs,NULL,&xsrefe,gpcoord[gpid]);
        double fac = gpweight * drs * curvefac;   // integration factor
        // distribute over element load vector
        for (int nodid=0; nodid < 2; ++nodid) {
          for(int dim=0; dim < NUMDIM_SOH8; ++dim) {
            elevec1[nodid*NUMDIM_SOH8 + dim] += funct[nodid] * (*onoff)[dim] * (*val)[dim] * fac;
          }
        }
      }
      break;
    default:
      dserror("Unknown type of LineNeumann load");
    break;
    }

  }
    return 0;
}

/*----------------------------------------------------------------------*
 * Evaluate sqrt of determinant of metric at gp (private)      maf 05/07*
 * ---------------------------------------------------------------------*/
void DRT::Elements::Soh8Line::soh8_line_integ(
      vector<double>* funct,                 // (o) shape functions
      double* sqrtdetg,                      // (o) pointer to sqrt of det(g)
      vector<double>* unrm,                  // (o) unit normal
      const Epetra_SerialDenseMatrix* xs,    // (i) element coords
      const double r)                        // (i) coord in r-direction
{
  // shape functions for 2 nodes
  (*funct)[0] = 0.5 * (1.0-r);
  (*funct)[1] = 0.5 * (1.0+r);
  // derivatives of 2 shape functions 
  Epetra_SerialDenseMatrix deriv(2,1);
  deriv(0,0) = -0.5;
  deriv(1,0) =  0.5;

  // compute dXYZ / drs
  Epetra_SerialDenseMatrix dxyzdrs(1,3);
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
  Epetra_SerialDenseMatrix metrictensor(1,1);
  metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,1.0);
  (*sqrtdetg) = sqrt( metrictensor(0,0));
  if (unrm != NULL){
    dserror(" unit normal on line not implemented!");
  }

  return;
}



#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOH8
