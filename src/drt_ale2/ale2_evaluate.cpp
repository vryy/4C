#ifdef D_ALE
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "ale2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"

using namespace DRT::Utils;

extern "C"
{
#include "../headers/standardtypes.h"
}

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::Elements::Ale2::Evaluate(ParameterList& params,
                                  DRT::Discretization&      discretization,
                                  vector<int>&              lm,
                                  Epetra_SerialDenseMatrix& elemat1,
                                  Epetra_SerialDenseMatrix& elemat2,
                                  Epetra_SerialDenseVector& elevec1,
                                  Epetra_SerialDenseVector& elevec2,
                                  Epetra_SerialDenseVector& elevec3)
{
  DRT::Elements::Ale2::ActionType act = Ale2::none;

  // get the action required
  string action = params.get<string>("action","none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_ale_lin_stiff")
    act = Ale2::calc_ale_lin_stiff;
  else
    dserror("Unknown type of action for Ale2");

  // get the material
  MATERIAL* actmat = &(mat[material_-1]);

  switch (act)
  {
  case calc_ale_lin_stiff:
  {
    //RefCountPtr<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
    //vector<double> my_dispnp(lm.size());
    //DRT::Utils::ExtractMyValues(*dispnp,my_dispnp,lm);

    static_ke(lm,&elemat1,&elevec1,actmat,params);

    break;
  }
  default:
    dserror("Unknown type of action for Ale2");
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      gammi 04/07|
 |                                                                      |
 |  The function is just a dummy. For the ale elements, the           |
 |  integration of the volume neumann loads takes place in the element. |
 |  We need it there for the stabilisation terms!                       |
 *----------------------------------------------------------------------*/
int DRT::Elements::Ale2::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


void DRT::Elements::Ale2::static_ke(vector<int>&              lm,
                                    Epetra_SerialDenseMatrix* sys_mat,
                                    Epetra_SerialDenseVector* residual,
                                    struct _MATERIAL*         material,
                                    ParameterList&            params)
{
  const int iel = NumNode();
  const int nd  = 2 * iel;
  const DiscretizationType distype = this->Shape();

  Epetra_SerialDenseMatrix xyze(2,iel);

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
  }

  /*----------------------------------------- declaration of variables ---*/
  Epetra_SerialDenseVector  funct(iel);
  Epetra_SerialDenseMatrix 	deriv(2,iel);
  Epetra_SerialDenseMatrix 	xjm(2,2);
  Epetra_SerialDenseMatrix 	xji(2,2);
  Epetra_SerialDenseMatrix 	bop(3,2*iel);
  Epetra_SerialDenseMatrix 	d(4,4);

  // gaussian points
  const GaussRule2D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints2D  intpoints = getIntegrationPoints2D(gaussrule);

  // integration loops
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
      const double e1 = intpoints.qxg[iquad][0];
      const double e2 = intpoints.qxg[iquad][1];

      // shape functions and their derivatives
      DRT::Utils::shape_function_2D(funct,e1,e2,distype);
      DRT::Utils::shape_function_2D_deriv1(deriv,e1,e2,distype);

      // compute jacobian matrix

      // determine jacobian at point r,s,t
      for (int i=0; i<2; i++)
      {
        for (int j=0; j<2; j++)
        {
          double dum=0.;
          for (int l=0; l<iel; l++)
          {
            dum += deriv(i,l)*xyze(j,l);
          }
          xjm(i,j)=dum;
        }
      }

      // determinant of jacobian
      const double det = xjm(0,0)*xjm(1,1) - xjm(0,1)*xjm(1,0);
      const double fac = intpoints.qwgt[iquad]*det;

      // calculate operator B

      // inverse of jacobian
      const double dum=1.0/det;
      xji(0,0) = xjm(1,1)* dum;
      xji(0,1) =-xjm(0,1)* dum;
      xji(1,0) =-xjm(1,0)* dum;
      xji(1,1) = xjm(0,0)* dum;

      // get operator b of global derivatives
      for (int i=0; i<iel; i++)
      {
        const int node_start = i*2;

        const double hr   = deriv(0,i);
        const double hs   = deriv(1,i);

        const double h1 = xji(0,0)*hr + xji(0,1)*hs;
        const double h2 = xji(1,0)*hr + xji(1,1)*hs;

        bop(0,node_start+0) = h1 ;
        bop(0,node_start+1) = 0.0;
        bop(1,node_start+0) = 0.0;
        bop(1,node_start+1) = h2 ;
        bop(2,node_start+0) = h2;
        bop(2,node_start+1) = h1;
      }

      // call material law

      const double ym  = material->m.stvenant->youngs;
      const double pv  = material->m.stvenant->possionratio;

      // plane strain, rotational symmetry
      const double c1=ym/(1.0+pv);
      const double b1=c1*pv/(1.0-2.0*pv);
      const double a1=b1+c1;

      d(0,0)=a1;
      d(0,1)=b1;
      d(0,2)=0.;
      d(0,3)=b1;

      d(1,0)=b1;
      d(1,1)=a1;
      d(1,2)=0.;
      d(1,3)=b1;

      d(2,0)=0.;
      d(2,1)=0.;
      d(2,2)=c1/2.;
      d(2,3)=0.;

      d(3,0)=b1;
      d(3,1)=b1;
      d(3,2)=0.;
      d(3,3)=a1;

      for (int j=0; j<nd; j++)
      {
        double db[3];
        for (int k=0; k<3; k++)
        {
          db[k] = 0.0;
          for (int l=0; l<3; l++)
          {
            db[k] += d(k,l)*bop(l,j)*fac ;
          }
        }
        for (int i=0; i<nd; i++)
        {
          double dum = 0.0;
          for (int m=0; m<3; m++)
          {
            dum = dum + bop(m,i)*db[m] ;
          }
          (*sys_mat)(i,j) += dum ;
        }
      }
  }
}


//=======================================================================
//=======================================================================

int DRT::Elements::Ale2Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif
#endif
#endif
