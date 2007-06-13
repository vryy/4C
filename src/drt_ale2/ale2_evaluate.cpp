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

  Epetra_SerialDenseMatrix xyze(2,iel);

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
  }

  /*----------------------------------------- declaration of variables ---*/
  vector<double> 		funct(iel);
  Epetra_SerialDenseMatrix 	deriv(2,iel);
  Epetra_SerialDenseMatrix 	xjm(2,2);
  Epetra_SerialDenseMatrix 	xji(2,2);
  Epetra_SerialDenseMatrix 	bop(3,2*iel);
  Epetra_SerialDenseMatrix 	d(4,4);
  ALE2_DATA                     data;

  double         		e1, e2;
  double         		facr=0.0, facs=0.0;
  int                           nir, nis;
  int                           intc;

  // gaussian points
  intg(data);

  // integration parameters
  switch (Shape())
  {
  case quad4: case quad8: case quad9:  /* --> quad - element */
    /* initialise integration */
    nir = ngp_[0];
    nis = ngp_[1];
    break;
  case tri6: /* --> tri - element */
  case tri3:
    /* initialise integration */
    nir  = ngp_[0];
    nis  = 1;
    intc = ngp_[1];
    break;
  default:
    dserror("typ unknown!");
  }

  // integration loops
  for (int lr=0; lr<nir; lr++)
  {
    for (int ls=0; ls<nis; ls++)
    {
      switch (Shape())
      {
      case quad4: case quad8: case quad9:   /* --> quad - element */
        e1   = data.qxg [lr][nir-1];
        facr = data.qwgt[lr][nir-1];
        e2   = data.qxg [ls][nis-1];
        facs = data.qwgt[ls][nis-1];
        break;
      case tri3: case tri6:   /* --> tri - element */
        e1   = data.txgr[lr][intc];
        facr = data.twgt[lr][intc];
        e2   = data.txgs[lr][intc];
        facs = 1.;
        break;
      default:
        facr = facs = 0.0;
        e1 = e2 = 0.0;
        dserror("typ unknown!");
      }

      // shape functions and their derivatives
      funct_deriv(funct,deriv,e1,e2,Shape(),1);

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
      double det = xjm(0,0)*xjm(1,1) - xjm(0,1)*xjm(1,0);

      double fac = facr * facs * det;

      // calculate operator B

      // inverse of jacobian
      double dum=1.0/det;
      xji(0,0) = xjm(1,1)* dum;
      xji(0,1) =-xjm(0,1)* dum;
      xji(1,0) =-xjm(1,0)* dum;
      xji(1,1) = xjm(0,0)* dum;

      // get operator b of global derivatives
      for (int i=0; i<iel; i++)
      {
        int node_start = i*2;

        double hr   = deriv(0,i);
        double hs   = deriv(1,i);

        double h1 = xji(0,0)*hr + xji(0,1)*hs;
        double h2 = xji(1,0)*hr + xji(1,1)*hs;

        bop(0,node_start+0) = h1 ;
        bop(0,node_start+1) = 0.0;
        bop(1,node_start+0) = 0.0;
        bop(1,node_start+1) = h2 ;
        bop(2,node_start+0) = h2;
        bop(2,node_start+1) = h1;
      }

      // call material law

      double ym  = material->m.stvenant->youngs;
      double pv  = material->m.stvenant->possionratio;

      // plane strain, rotational symmetry
      double c1=ym/(1.0+pv);
      double b1=c1*pv/(1.0-2.0*pv);
      double a1=b1+c1;

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
}



void DRT::Elements::Ale2::intg(ALE2_DATA& data)
{
  for (int i=0; i<13; i++) /* loop over all integration points    */
  {
    for (int k=0; k<11; k++) /* loop integration cases          */
    {
      /* set coordinates (r,s) coordinates of integration point     */
      data.txgr[i][k] = 0.;
      data.txgs[i][k] = 0.;

      /* innitialise the vector of gaussweights */
      data.twgt[i][k] = 0.;
    }
  }
  /* inintialise quadX arrays */
  for (int i=0; i<6; i++) /* loop over integration points       */
  {
    for (int k=0; k<6; k++) /* loop integration cases         */
    {
      /* set one coordinate of integration points --- the rest is
       * 'symmetric'                                               */
      data.qxg [i][k] = 0.;

      /* innitialise the vector of gaussweights */
      data.qwgt[i][k] = 0.;
    }
  }

/*----------------------------------------------------------------------*
  |     INTEGRATION PARAMETERS FOR    R E C T A N G U L A R   ELEMENTS   |
  |     GAUSS SAMPLING POINTS  AT     R/S-COORDINATES     RESPECTIVELY   |
  |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
  |    data.qxg[i][j]                                                    |
  |   data.qwgt[i][j]:  i+1 - actual number of gausspoint                |
  |                     j+1 - total number of gausspoints                |
  *----------------------------------------------------------------------*/
/* coordinates for two gauss points */
  data.qxg[0][1]  =  -0.5773502691896;
  data.qxg[1][1]  =  -data.qxg[0][1] ;
/* coordinates for three gauss points */
  data.qxg[0][2]  =  -0.7745966692415;
  data.qxg[2][2]  =  -data.qxg[0][2] ;
/* coordinates for four gauss points */
  data.qxg[0][3]  =  -0.8611363115941;
  data.qxg[1][3]  =  -0.3399810435849;
  data.qxg[2][3]  =  -data.qxg[1][3] ;
  data.qxg[3][3]  =  -data.qxg[0][3] ;
/* coordinates for five gauss points */
  data.qxg[0][4]  =  -0.9061798459387;
  data.qxg[1][4]  =  -0.5384693101057;
  data.qxg[3][4]  =  -data.qxg[1][4] ;
  data.qxg[4][4]  =  -data.qxg[0][4] ;
/* coordinates for six gauss points */
  data.qxg[0][5]  =  -0.9324695142032;
  data.qxg[1][5]  =  -0.6612093864663;
  data.qxg[2][5]  =  -0.2386191860832;
  data.qxg[3][5]  =  -data.qxg[2][5] ;
  data.qxg[4][5]  =  -data.qxg[1][5] ;
  data.qxg[5][5]  =  -data.qxg[0][5] ;

/* weights for one gauss points */
  data.qwgt[0][0] =  TWO             ;
/* weights for two gauss points */
  data.qwgt[0][1] =  ONE             ;
  data.qwgt[1][1] =  ONE             ;
/* weights for three gauss points */
  data.qwgt[0][2] =  0.5555555555556 ;
  data.qwgt[1][2] =  0.8888888888889 ;
  data.qwgt[2][2] =  data.qwgt[0][2] ;
/* weights for four gauss points */
  data.qwgt[0][3] =  0.3478548451375 ;
  data.qwgt[1][3] =  0.6521451548625 ;
  data.qwgt[2][3] =  data.qwgt[1][3] ;
  data.qwgt[3][3] =  data.qwgt[0][3] ;
/* weights for five gauss points */
  data.qwgt[0][4] =  0.2369268850562 ;
  data.qwgt[1][4] =  0.4786286704994 ;
  data.qwgt[2][4] =  0.5688888888889 ;
  data.qwgt[3][4] =  data.qwgt[1][4] ;
  data.qwgt[4][4] =  data.qwgt[0][4] ;
/* weights for six gauss points */
  data.qwgt[0][5] =  0.1713244923792 ;
  data.qwgt[1][5] =  0.3607615730481 ;
  data.qwgt[2][5] =  0.4679139345727 ;
  data.qwgt[3][5] =  data.qwgt[2][5] ;
  data.qwgt[4][5] =  data.qwgt[1][5] ;
  data.qwgt[5][5] =  data.qwgt[0][5] ;


/*----------------------------------------------------------------------*
  |     INTEGRATION PARAMETERS FOR    T R I A N G U L A R     ELEMENTS   |
  |     GAUSS SAMPLING POINTS  AT     R/S-COORDINATES     RESPECTIVELY   |
  |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
  |   data.txgr[i][j]                                                    |
  |  data.twgts[i][j]:  i+1 - actual number of gausspoint                |
  |                     j+1 - number for integration case (from input)   |
  *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
  |    GAUSS INTEGRATION         1 SAMPLING POINT, DEG.OF PRECISION 1    |
  |                              CASE 0                                  |
  *----------------------------------------------------------------------*/
  data.txgr[0][0]    =  1./3. ;
  data.txgs[0][0]    =  1./3. ;

  data.twgt[0][0]    =  1./2. ;
/*----------------------------------------------------------------------*
  |    GAUSS INTEGRATION        3 SAMPLING POINTS, DEG.OF PRECISION 2    |
  |                             CASE 1                                   |
  *----------------------------------------------------------------------*/
  data.txgr[0][1]    =  1./2.  ;
  data.txgr[1][1]    =  1./2.  ;
  data.txgr[2][1]    =  0. ;
  data.txgs[0][1]    =  0. ;
  data.txgs[1][1]    =  1./2.  ;
  data.txgs[2][1]    =  1./2.  ;

  data.twgt[0][1]    =  1./6.  ;
  data.twgt[1][1]    =  1./6.  ;
  data.twgt[2][1]    =  1./6.  ;
/*----------------------------------------------------------------------*
  |    ALT.GAUSS INTEGRATION    3 SAMPLING POINTS, DEG.OF PRECISION 2    |
  |                             CASE 2                                   |
  *----------------------------------------------------------------------*/
  data.txgr[0][2]    =  1./6.  ;
  data.txgr[1][2]    =  2./3.  ;
  data.txgr[2][2]    =  1./6.  ;
  data.txgs[0][2]    =  1./6.  ;
  data.txgs[1][2]    =  1./6.  ;
  data.txgs[2][2]    =  2./3.  ;

  data.twgt[0][2]    =  1./6.  ;
  data.twgt[1][2]    =  1./6.  ;
  data.twgt[2][2]    =  1./6.  ;
/*----------------------------------------------------------------------*
  |    GAUSS INTEGRATION        4 SAMPLING POINTS, DEG.OF PRECISION 3    |
  |                             CASE 3                                   |
  *----------------------------------------------------------------------*/
  data.txgr[0][3]    =  0.2                ;
  data.txgr[1][3]    =  0.6                ;
  data.txgr[2][3]    =  0.2                ;
  data.txgr[3][3]    =  1./3.                ;
  data.txgs[0][3]    =  0.2                ;
  data.txgs[1][3]    =  0.2                ;
  data.txgs[2][3]    =  0.6                ;
  data.txgs[3][3]    =  1./3.                ;

  data.twgt[0][3]    =  0.2604166666667    ;
  data.twgt[1][3]    =  data.twgt[0][2]    ;
  data.twgt[2][3]    =  data.twgt[0][2]    ;
  data.twgt[3][3]    = -0.28125            ;
/*----------------------------------------------------------------------*
  |    GAUSS INTEGRATION        6 SAMPLING POINTS, DEG.OF PRECISION 4    |
  |                             CASE 4                                   |
  *----------------------------------------------------------------------*/
  data.txgr[0][4]    =  0.0915762135098	;
  data.txgr[1][4]    =  0.8168475729805	;
  data.txgr[2][4]    =  0.0915762135098 	;
  data.txgr[3][4]    =  0.4459484909160	;
  data.txgr[4][4]    =  0.4459484909160 	;
  data.txgr[5][4]    =  0.1081030181681	;
  data.txgs[0][4]    =  0.0915762135098 	;
  data.txgs[1][4]    =  0.0915762135098 	;
  data.txgs[2][4]    =  0.8168475729805 	;
  data.txgs[3][4]    =  0.1081030181681 	;
  data.txgs[4][4]    =  0.4459484909160 	;
  data.txgs[5][4]    =  0.4459484909160 	;

  data.twgt[0][4]   =  0.0549758718277	;
  data.twgt[1][4]   =  0.0549758718277	;
  data.twgt[2][4]   =  0.0549758718277	;
  data.twgt[3][4]   =  0.1116907948390	;
  data.twgt[4][4]   =  0.1116907948390	;
  data.twgt[5][4]   =  0.1116907948390	;
/*----------------------------------------------------------------------*
  |    ALT.GAUSS INTEGRATION    6 SAMPLING POINTS, DEG.OF PRECISION 3    |
  |                             CASE 5                                   |
  *----------------------------------------------------------------------*/
  data.txgr[0][5]    =  0.1090390090729	;
  data.txgr[1][5]    =  0.2319333685530	;
  data.txgr[2][5]    =  0.6590276223741	;
  data.txgr[3][5]    =  0.0915762135098 	;
  data.txgr[4][5]    =  0.8168475729805 	;
  data.txgr[5][5]    =  0.0915762135098 	;
  data.txgs[0][5]    =  0.8168475729805 	;
  data.txgs[1][5]    =  0.0915762135098 	;
  data.txgs[2][5]    =  0.0915762135098 	;
  data.txgs[3][5]    =  0.8168475729805 	;
  data.txgs[4][5]    =  0.0915762135098 	;
  data.txgs[5][5]    =  0.0915762135098 	;



  data.twgt[0][5]   =  0.0833333333333	;
  data.twgt[1][5]   =  0.0833333333333	;
  data.twgt[2][5]   =  0.0833333333333	;
  data.twgt[3][5]   =  0.0833333333333	;
  data.twgt[4][5]   =  0.0833333333333	;
  data.twgt[5][5]   =  0.0833333333333	;
/*----------------------------------------------------------------------*
  |    GAUSS INTEGRATION        7 SAMPLING POINTS, DEG.OF PRECISION 5    |
  |                             CASE 6                                   |
  *----------------------------------------------------------------------*/
  data.txgr[0][6]    =  0.1012865073235 ;
  data.txgr[1][6]    =  0.4701420641051 ;
  data.txgr[2][6]    =  0.7974269853531 ;
  data.txgr[3][6]    =  data.txgr[1][4]       ;
  data.txgr[4][6]    =  data.txgr[0][4]       ;
  data.txgr[5][6]    =  0.0597158717898 ;
  data.txgr[6][6]    =  1./3.	      ;
  data.txgs[0][6]    =  data.txgr[0][4]       ;
  data.txgs[1][6]    =  data.txgr[5][4]       ;
  data.txgs[2][6]    =  data.txgr[0][4]       ;
  data.txgs[3][6]    =  data.txgr[1][4]       ;
  data.txgs[4][6]    =  data.txgr[2][4]       ;
  data.txgs[5][6]    =  data.txgr[1][4]       ;
  data.txgs[6][6]    =  1./3.	      ;

  data.twgt[0][6]    =  0.0629695902724 ;
  data.twgt[1][6]    =  0.0661970763943 ;
  data.twgt[2][6]    =  data.twgt[0][4]      ;
  data.twgt[3][6]    =  data.twgt[1][4]      ;
  data.twgt[4][6]    =  data.twgt[0][4]      ;
  data.twgt[5][6]    =  data.twgt[1][4]      ;
  data.twgt[6][6]    =  0.1125	      ;
/*----------------------------------------------------------------------*
  |    ALT.GAUSS INTEGRATION    7 SAMPLING POINTS, DEG.OF PRECISION 4    |
  |                             CASE 7                                   |
  *----------------------------------------------------------------------*/
  data.txgr[0][7]    =  0.2379323664724 ;
  data.txgr[1][7]    =  0.7367124989684 ;
  data.txgr[2][7]    =  data.txgr[1][4] ;
  data.txgr[3][7]    =  data.txgr[0][4] ;
  data.txgr[4][7]    =  0.0253551345591 ;
  data.txgr[5][7]    =  data.txgr[4][4] ;
  data.txgr[6][7]    =  1./3.	            ;
  data.txgs[0][7]    =  data.txgr[4][4] ;
  data.txgs[1][7]    =  data.txgr[4][4] ;
  data.txgs[2][7]    =  data.txgr[0][4] ;
  data.txgs[3][7]    =  data.txgr[1][4] ;
  data.txgs[4][7]    =  data.txgr[1][4] ;
  data.txgs[5][7]    =  data.txgr[0][4] ;
  data.txgs[6][7]    =  1./3.	            ;

  data.twgt[0][7]    =  0.0520833333333 ;
  data.twgt[1][7]    =  data.twgt[0][4] ;
  data.twgt[2][7]    =  data.twgt[0][4] ;
  data.twgt[3][7]    =  data.twgt[0][4] ;
  data.twgt[4][7]    =  data.twgt[0][4] ;
  data.twgt[5][7]    =  data.twgt[0][4] ;
  data.twgt[6][7]    =  0.1875	    ;
/*----------------------------------------------------------------------*
  |    GAUSS INTEGRATION        9 SAMPLING POINTS, DEG.OF PRECISION 5    |
  |                             CASE 8                                   |
  *----------------------------------------------------------------------*/
  data.txgr[0][8]    =  0.1654099273898     ;
  data.txgr[1][8]    =  0.4375252483834     ;
  data.txgr[2][8]    =  0.7971126518601     ;
  data.txgr[3][8]    =  data.txgr[2][5]     ;
  data.txgr[4][8]    =  data.txgr[1][5]     ;
  data.txgr[5][8]    =  data.txgr[0][5]     ;
  data.txgr[6][8]    =  0.0374774207501     ;
  data.txgr[7][8]    =  0.1249495032332     ;
  data.txgr[8][8]    =  data.txgr[6][5]     ;
  data.txgs[0][8]    =  data.txgr[6][5]     ;
  data.txgs[1][8]    =  data.txgr[7][5]     ;
  data.txgs[2][8]    =  data.txgr[6][5]     ;
  data.txgs[3][8]    =  data.txgr[0][5]     ;
  data.txgs[4][8]    =  data.txgr[1][5]     ;
  data.txgs[5][8]    =  data.txgr[2][5]     ;
  data.txgs[6][8]    =  data.txgr[2][5]     ;
  data.txgs[7][8]    =  data.txgr[1][5]     ;
  data.txgs[8][8]    =  data.txgr[0][5]     ;

  data.twgt[0][8]    =  0.0318457071431     ;
  data.twgt[1][8]    =  0.1029752523804     ;
  data.twgt[2][8]    =  data.twgt[0][5]     ;
  data.twgt[3][8]    =  data.twgt[0][5]     ;
  data.twgt[4][8]    =  data.twgt[1][5]     ;
  data.twgt[5][8]    =  data.twgt[0][5]     ;
  data.twgt[6][8]    =  data.twgt[0][5]     ;
  data.twgt[7][8]    =  data.twgt[1][5]     ;
  data.twgt[8][8]    =  data.twgt[0][5]     ;
  /*----------------------------------------------------------------------*
    |    GAUSS INTEGRATION       12 SAMPLING POINTS, DEG.OF PRECISION 6    |
    |                            CASE 9                                    |
    *----------------------------------------------------------------------*/
  data.txgr[ 0][9]   =  0.0630890144915     ;
  data.txgr[ 1][9]   =  0.3103524510338     ;
  data.txgr[ 2][9]   =  0.6365024991214     ;
  data.txgr[ 3][9]   =  0.8738219710170     ;
  data.txgr[ 4][9]   =  data.txgr[ 2][6]    ;
  data.txgr[ 5][9]   =  data.txgr[ 1][6]    ;
  data.txgr[ 6][9]   =  data.txgr[ 0][6]    ;
  data.txgr[ 7][9]   =  0.0531450498448     ;
  data.txgr[ 8][9]   =  data.txgr[ 7][6]    ;
  data.txgr[ 9][9]   =  0.2492867451709     ;
  data.txgr[10][9]   =  0.5014265096582     ;
  data.txgr[11][9]   =  data.txgr[ 9][6]    ;
  data.txgs[ 0][9]   =  data.txgr[ 0][6]    ;
  data.txgs[ 1][9]   =  data.txgr[ 7][6]    ;
  data.txgs[ 2][9]   =  data.txgr[ 7][6]    ;
  data.txgs[ 3][9]   =  data.txgr[ 0][6]    ;
  data.txgs[ 4][9]   =  data.txgr[ 1][6]    ;
  data.txgs[ 5][9]   =  data.txgr[ 2][6]    ;
  data.txgs[ 6][9]   =  data.txgr[ 3][6]    ;
  data.txgs[ 7][9]   =  data.txgr[ 2][6]    ;
  data.txgs[ 8][9]   =  data.txgr[ 1][6]    ;
  data.txgs[ 9][9]   =  data.txgr[ 9][6]    ;
  data.txgs[10][9]   =  data.txgr[ 9][6]    ;
  data.txgs[11][9]   =  data.txgr[10][6]    ;
  data.twgt[ 0][9]   =  0.0254224531851     ;
  data.twgt[ 1][9]   =  0.0414255378092     ;
  data.twgt[ 2][9]   =  data.twgt[ 1][6]    ;
  data.twgt[ 3][9]   =  data.twgt[ 0][6]    ;
  data.twgt[ 4][9]   =  data.twgt[ 1][6]    ;
  data.twgt[ 5][9]   =  data.twgt[ 1][6]    ;
  data.twgt[ 6][9]   =  data.twgt[ 0][6]    ;
  data.twgt[ 7][9]   =  data.twgt[ 1][6]    ;
  data.twgt[ 8][9]   =  data.twgt[ 1][6]    ;
  data.twgt[ 9][9]   =  0.0583931378632     ;
  data.twgt[10][9]   =  data.twgt[ 9][6]    ;
  data.twgt[11][6]   =  data.twgt[ 9][6]    ;
/*----------------------------------------------------------------------*
  |    GAUSS INTEGRATION       13 SAMPLING POINTS, DEG.OF PRECISION 7    |
  |                            CASE 10                                   |
  *----------------------------------------------------------------------*/
  data.txgr[ 0][10]  =  0.0651301029022     ;
  data.txgr[ 1][10]  =  0.3128654960049     ;
  data.txgr[ 2][10]  =  0.6384441885698     ;
  data.txgr[ 3][10]  =  0.8697397941956     ;
  data.txgr[ 4][10]  =  data.txgr[ 2][7]    ;
  data.txgr[ 5][10]  =  data.txgr[ 1][7]    ;
  data.txgr[ 6][10]  =  data.txgr[ 0][7]    ;
  data.txgr[ 7][10]  =  0.0486903154253     ;
  data.txgr[ 8][10]  =  data.txgr[ 7][7]    ;
  data.txgr[ 9][10]  =  0.2603459660790     ;
  data.txgr[10][10]  =  0.4793080678419     ;
  data.txgr[11][10]  =  data.txgr[ 9][7]    ;
  data.txgr[12][10]  =  1./3.                 ;
  data.txgs[ 0][10]  =  data.txgr[ 0][7]    ;
  data.txgs[ 1][10]  =  data.txgr[ 7][7]    ;
  data.txgs[ 2][10]  =  data.txgr[ 7][7]    ;
  data.txgs[ 3][10]  =  data.txgr[ 0][7]    ;
  data.txgs[ 4][10]  =  data.txgr[ 1][7]    ;
  data.txgs[ 5][10]  =  data.txgr[ 2][7]    ;
  data.txgs[ 6][10]  =  data.txgr[ 3][7]    ;
  data.txgs[ 7][10]  =  data.txgr[ 2][7]    ;
  data.txgs[ 8][10]  =  data.txgr[ 1][7]    ;
  data.txgs[ 9][10]  =  data.txgr[ 9][7]    ;
  data.txgs[10][10]  =  data.txgr[ 9][7]    ;
  data.txgs[11][10]  =  data.txgr[10][7]    ;
  data.txgs[12][10]  =  1./3.                 ;
  data.twgt[ 0][10]  =  0.0266736178044     ;
  data.twgt[ 1][10]  =  0.0385568804451     ;
  data.twgt[ 2][10]  =  data.twgt[ 1][7]    ;
  data.twgt[ 3][10]  =  data.twgt[ 0][7]    ;
  data.twgt[ 4][10]  =  data.twgt[ 1][7]    ;
  data.twgt[ 5][10]  =  data.twgt[ 1][7]    ;
  data.twgt[ 6][10]  =  data.twgt[ 0][7]    ;
  data.twgt[ 7][10]  =  data.twgt[ 1][7]    ;
  data.twgt[ 8][10]  =  data.twgt[ 1][7]    ;
  data.twgt[ 9][10]  =  0.0878076287166     ;
  data.twgt[10][10]  =  data.twgt[ 9][7]    ;
  data.twgt[11][10]  =  data.twgt[ 9][7]    ;
  data.twgt[12][10]  = -0.0747850222338     ;
}


void DRT::Elements::Ale2::funct_deriv(
    vector<double>& funct,
    Epetra_SerialDenseMatrix& deriv,
    double      r,
    double      s,
    DiscretizationType     typ,
    int         option
    )
{
  const DOUBLE   q18 = 1.0/8.0;
  const DOUBLE   q14 = 1.0/4.0;
  DOUBLE         rp,sp,rm,sm,rrm,ssm;
  DOUBLE         t1,t2,t3,t4;

  /* if option ==0 only funtion evaluation, if option==1 also derivatives */
  rp  = 1.0+r;
  rm  = 1.0-r;
  sp  = 1.0+s;
  sm  = 1.0-s;
  rrm = 1.0-r*r;
  ssm = 1.0-s*s;

  switch (typ)
  {
  case quad4: /* LINEAR shape functions for quad4 and their natural
           *                                          derivatives ----*/
  {
    /*--------------------------------------------- form basic values */
    double rp=1.+r;
    double rm=1.-r;
    double sp=1.+s;
    double sm=1.-s;

    funct[0]=q14*rp*sp;
    funct[1]=q14*rm*sp;
    funct[2]=q14*rm*sm;
    funct[3]=q14*rp*sm;

    if (option==1)
    {
      deriv(0,0)= q14*sp;
      deriv(1,0)= q14*rp;

      deriv(0,1)=-q14*sp;
      deriv(1,1)= q14*rm;

      deriv(0,2)=-q14*sm;
      deriv(1,2)=-q14*rm;

      deriv(0,3)= q14*sm;
      deriv(1,3)=-q14*rp;
    }
    break;
  }
  case quad8: /* QUADRATIC shape functions for quadrilaterals without
             central node and their natural derivatives (serendipity) */
  {
    double rp=1.+r;
    double rm=1.-r;
    double sp=1.+s;
    double sm=1.-s;
    double r2=1.-r*r;
    double s2=1.-s*s;

    funct[4]=1./2.*r2*sp;
    funct[5]=1./2.*rm*s2;
    funct[6]=1./2.*r2*sm;
    funct[7]=1./2.*rp*s2;
    funct[0]=q14*rp*sp-1./2.*(funct[4]+funct[7]);
    funct[1]=q14*rm*sp-1./2.*(funct[4]+funct[5]);
    funct[2]=q14*rm*sm-1./2.*(funct[5]+funct[6]);
    funct[3]=q14*rp*sm-1./2.*(funct[6]+funct[7]);

    if (option==1)
    {
      deriv(0,0)= q14*sp;
      deriv(1,0)= q14*rp;

      deriv(0,1)=-q14*sp;
      deriv(1,1)= q14*rm;

      deriv(0,2)=-q14*sm;
      deriv(1,2)=-q14*rm;

      deriv(0,3)= q14*sm;
      deriv(1,3)=-q14*rp;

      deriv(0,4)=-1.*r*sp;
      deriv(1,4)= 1./2.*r2;

      deriv(0,5)=-1./2.*  s2;
      deriv(1,5)=-1.*rm*s;

      deriv(0,6)=-1.*r*sm;
      deriv(1,6)=-1./2.*r2;

      deriv(0,7)= 1./2.*  s2;
      deriv(1,7)=-1.*rp*s;

      deriv(0,0)-= 1./2.*(deriv(0,4)+deriv(0,7));
      deriv(1,0)-= 1./2.*(deriv(1,4)+deriv(1,7));

      for(int i=1;i<4;i++)
      {
        int ii=i+3;
        deriv(0,i) -= 1./2.*(deriv(0,ii)+deriv(0,ii+1));
        deriv(1,i) -= 1./2.*(deriv(1,ii)+deriv(1,ii+1));
      }
    }
    break;
  }
  case quad9: /* full QUADRATIC shape functions for quadrilaterals with
             central node and their natural derivatives */
  {
/*--------------------------------------------------- form basic values */
    double rp=1.+r;
    double rm=1.-r;
    double sp=1.+s;
    double sm=1.-s;
    double r2=1.-r*r;
    double s2=1.-s*s;
    double rh=1./2.*r;
    double sh=1./2.*s;
    double rs=rh*sh;
    double rhp=r+1./2.;
    double rhm=r-1./2.;
    double shp=s+1./2.;
    double shm=s-1./2.;

    funct[0]= rs*rp*sp;
    funct[1]=-rs*rm*sp;
    funct[2]= rs*rm*sm;
    funct[3]=-rs*rp*sm;
    funct[4]= sh*sp*r2;
    funct[5]=-rh*rm*s2;
    funct[6]=-sh*sm*r2;
    funct[7]= rh*rp*s2;
    funct[8]= r2*s2;

    if (option==1)
    {
      deriv(0,0)= rhp*sh*sp;
      deriv(1,0)= shp*rh*rp;

      deriv(0,1)= rhm*sh*sp;
      deriv(1,1)=-shp*rh*rm;

      deriv(0,2)=-rhm*sh*sm;
      deriv(1,2)=-shm*rh*rm;

      deriv(0,3)=-rhp*sh*sm;
      deriv(1,3)= shm*rh*rp;

      deriv(0,4)=-2.*r*sh*sp;
      deriv(1,4)= shp*r2;

      deriv(0,5)= rhm*s2;
      deriv(1,5)= 2.*s*rh*rm;

      deriv(0,6)= 2.*r*sh*sm;
      deriv(1,6)= shm*r2;

      deriv(0,7)= rhp*s2;
      deriv(1,7)=-2.*s*rh*rp;

      deriv(0,8)=-2.*r*s2;
      deriv(1,8)=-2.*s*r2;
    } /* endif (icode>1) */
    break;
  }
  case tri3: /* LINEAR shape functions for triangles and their natural
           *                                         derivatives -----*/
  {
    /*------------------------------------------- form basic values */
    funct[0]=1.-r-s;
    funct[1]=r;
    funct[2]=s;

    if (option==1)
    {
      deriv(0,0)=-1.;
      deriv(1,0)=-1.;
      deriv(0,1)= 1.;
      deriv(1,1)=0.;
      deriv(0,2)=0.;
      deriv(1,2)= 1.;
    } /* endif (icode>1) */
    break;
  }
  case tri6: /* QUADRATIC shape functions for triangles and their natural
           *                                             derivatives -*/
  {
    /*------------------------------------------- form basic values */
    double rr=r*r;
    double ss=s*s;
    double rs=r*s;

    funct[0]=(1.-2.*r-2.*s)*(1.-r-s);
    funct[1]=2.*rr-r;
    funct[2]=2.*ss-s;
    funct[3]=4.*(r-rr-rs);
    funct[4]=4.*rs;
    funct[5]=4.*(s-rs-ss);

    if (option==1)
    {
      deriv(0,0)=-THREE+4.*(r+s);
      deriv(1,0)= deriv(0,0);

      deriv(0,1)= 4.*r-1.;
      deriv(1,1)= 0.;

      deriv(0,2)= 0.;
      deriv(1,2)= 4.*s-1.;

      deriv(0,3)= 4.*(1.-2.*r-s);
      deriv(1,3)=-4.*r;

      deriv(0,4)= 4.*s;
      deriv(1,4)= 4.*r;

      deriv(0,5)=-4.*s;
      deriv(1,5)= 4.*(1.-r-2.*s);
    } /* endif (icode>1) */
    break;
  }
  /*------------------------------------------------------------------*/
  default:
    dserror("distyp unknown\n");
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
