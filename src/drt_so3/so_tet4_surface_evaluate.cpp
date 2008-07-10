/*!----------------------------------------------------------------------*###
\file so_tet4_surface_evaluate.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
written by : Alexander Volf
			alexander.volf@mytum.de  
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_tet4.H"
#include "so_integrator.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"

/*----------------------------------------------------------------------*
 * Integrate a Surface Neumann boundary condition (public)     vlf 08/07*
 * ---------------------------------------------------------------------*/
int DRT::ELEMENTS::Sotet4Surface::EvaluateNeumann(ParameterList&           params,
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
  else if (*type == "neum_orthopressure") ltype = neum_orthopressure;
  else dserror("Unknown type of SurfaceNeumann condition");
  
  // get values and switches from the condition
  static const DRT::ELEMENTS::Integrator_tri3_1point tri3_int;
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

  // element geometry
  const int numnod = 3;
  Epetra_SerialDenseMatrix xsrefe(numnod,NUMDIM_SOTET4);  // material coord. of element
  for (int i=0; i<numnod; i++){
    xsrefe(i,0) = Nodes()[i]->X()[0];
    xsrefe(i,1) = Nodes()[i]->X()[1];
    xsrefe(i,2) = Nodes()[i]->X()[2];
  }

  Epetra_SerialDenseVector A(NUMDIM_SOTET4);
  Epetra_SerialDenseVector B(NUMDIM_SOTET4);
  Epetra_SerialDenseVector C(NUMDIM_SOTET4);
  
  /*
   * to compute the Jacobian we compute the area of the triangle 
   * with boundary vectors A,B of the triangle the area is |A x B| 
   */
  
  A(0)=xsrefe(1,0)-xsrefe(0,0);
  A(1)=xsrefe(1,1)-xsrefe(0,1);
  A(2)=xsrefe(1,2)-xsrefe(0,2);

  
  B(0)=xsrefe(2,0)-xsrefe(0,0);
  B(1)=xsrefe(2,1)-xsrefe(0,1);
  B(2)=xsrefe(2,2)-xsrefe(0,2);
  
  /* C = A x B  */
  C(0)=A(0)*B(1) - A(1)*B(0);
  C(1)=A(1)*B(2) - A(2)*B(1);
  C(2)=A(2)*B(0) - A(0)*B(2);
  
  /* detJ = |A x B|*/
  double detJ= C.Norm2()/2;

  /*
  ** Here, we integrate a 3-node surface with 1 Gauss Point
  */
  switch(ltype){
    case neum_live:{            // uniform load on reference configuration
      double fac = tri3_int.weights[0] * detJ * curvefac;   // integration factor
      //dserror("3 gauss points times 0.5 * g_1 \times g_2 correct?");
      // gauss parameters
      for (int gpid = 0; gpid < NUMGPT_SOTET4_FACE; gpid++) {    // loop over intergration points
        // get shape functions and derivatives of element surface
        // distribute over element load vector
        for (int nodid=0; nodid < NUMNOD_SOTET4_FACE; nodid++) {
          for(int dim=0; dim < NUMDIM_SOTET4; dim++) {
            elevec1[nodid*NUMDIM_SOTET4 + dim] +=
            	tri3_int.shapefct_gp[gpid](nodid) * (*onoff)[dim] * (*val)[dim] * fac;
          }
        }
      }
    }
    break;
    case neum_orthopressure:{   // orthogonal pressure on deformed config.
      dserror("orthopressure not implemented for tet4");
    }
    break;
    default:
      dserror("Unknown type of SurfaceNeumann load");
    break;
  }
  
  return 0;
} //Sotet4Surface::EvaluateNeumann(..)

#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOLID3
