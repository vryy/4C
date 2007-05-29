#ifdef D_ALE
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "ale3.H"
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
int DRT::Elements::Ale3::Evaluate(ParameterList& params,
                                  DRT::Discretization&      discretization,
                                  vector<int>&              lm,
                                  Epetra_SerialDenseMatrix& elemat1,
                                  Epetra_SerialDenseMatrix& elemat2,
                                  Epetra_SerialDenseVector& elevec1,
                                  Epetra_SerialDenseVector& elevec2,
                                  Epetra_SerialDenseVector& elevec3)
{
  DRT::Elements::Ale3::ActionType act = Ale3::none;

  // get the action required
  string action = params.get<string>("action","none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_ale_lin_stiff")
    act = Ale3::calc_ale_lin_stiff;
  else
    dserror("Unknown type of action for Ale3");

  // get the material
  MATERIAL* actmat = &(mat[material_-1]);

  switch(act)
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
    dserror("Unknown type of action for Ale3");
  }

  return 0;
}


extern "C"
{
  void dyn_facfromcurve(int actcurve,double T,double *fac);
}


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      gammi 04/07|
 |                                                                      |
 |  The function is just a dummy. For the ale elements, the           |
 |  integration of the volume neumann loads takes place in the element. |
 |  We need it there for the stabilisation terms!                       |
 *----------------------------------------------------------------------*/
int DRT::Elements::Ale3::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  return 0;
}


void DRT::Elements::Ale3::static_ke(vector<int>&              lm,
                                    Epetra_SerialDenseMatrix* sys_mat,
                                    Epetra_SerialDenseVector* residual,
                                    struct _MATERIAL*         material,
                                    ParameterList&            params)
{
  const int iel = NumNode();
  const int nd  = 3 * iel;

  Epetra_SerialDenseMatrix xyze(3,iel);

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
    xyze(2,i)=Nodes()[i]->X()[2];
  }

  /*----------------------------------------- declaration of variables ---*/
  vector<double> 		funct(iel);
  Epetra_SerialDenseMatrix 	deriv(3,iel);
  Epetra_SerialDenseMatrix 	xjm(3,3);
  Epetra_SerialDenseMatrix 	bop(6,3*iel);
  Epetra_SerialDenseMatrix 	D(6,6);
  ALE3_DATA                     data;

  double         		e1, e2, e3;
  double         		facr=0.0, facs=0.0, fact=0.0;
  int                           nir, nis, nit;
  double                        vol=0.;

  // gaussian points
  intg(data);

  // integration parameters
  switch (Shape())
  {
    case hex8:
    case hex20:
    case hex27:
      nir     = ngp_[0];
      nis     = ngp_[1];
      nit     = ngp_[2];
      break;
    case tet4:
    case tet10:
      nir = ngp_[0];
      nis = 1;
      nit = 1;
      break;
    default:
      dserror("unknown number of gaussian points in ale2_intg");
      break;
  }

  // integration loops
  for (int lr=0; lr<nir; lr++)
  {
    // gaussian point and weight at it
    e1   = data.xgpr[lr];
    facr = data.wgtr[lr];
    for (int ls=0; ls<nis; ls++)
    {
      // gaussian point and weight at it
      for (int lt=0; lt<nit; lt++)
      {
        // gaussian point and weight at it
        switch (Shape())
        {
          case hex8:
          case hex20:
          case hex27:
            e2   = data.xgps[ls];
            facs = data.wgts[ls];
            e3   = data.xgpt[lt];
            fact = data.wgtt[lt];
            break;
          case tet4:
          case tet10:
            e2   = data.xgps[lr];
            facs = 1.;
            e3   = data.xgpt[lr];
            fact = 1.;
            break;
          default:
            dserror("unknown number of gaussian points in ale2_intg");
            break;
        }

        // shape functions and their derivatives
        funct_deriv(funct,deriv,e1,e2,e3,Shape(),1);

        // compute jacobian matrix

        // determine jacobian at point r,s,t
        for (int i=0; i<3; i++)
        {
          for (int j=0; j<3; j++)
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
        double det = xjm(0,0)*xjm(1,1)*xjm(2,2) +
                     xjm(0,1)*xjm(1,2)*xjm(2,0) +
                     xjm(0,2)*xjm(1,0)*xjm(2,1) -
                     xjm(0,2)*xjm(1,1)*xjm(2,0) -
                     xjm(0,0)*xjm(1,2)*xjm(2,1) -
                     xjm(0,1)*xjm(1,0)*xjm(2,2);


        // calculate element volume
        vol += facr*facs*fact*det;

        double fac = facr * facs * fact * det;

        // calculate operator B

        // inverse of jacobian
        double x1r = xjm(0,0);
        double x2r = xjm(0,1);
        double x3r = xjm(0,2);
        double x1s = xjm(1,0);
        double x2s = xjm(1,1);
        double x3s = xjm(1,2);
        double x1t = xjm(2,0);
        double x2t = xjm(2,1);
        double x3t = xjm(2,2);

        double dum=1.0/det;

        double xi11=dum*(x2s*x3t - x2t*x3s);
        double xi12=dum*(x3r*x2t - x2r*x3t);
        double xi13=dum*(x2r*x3s - x3r*x2s);
        double xi21=dum*(x3s*x1t - x3t*x1s);
        double xi22=dum*(x1r*x3t - x3r*x1t);
        double xi23=dum*(x3r*x1s - x1r*x3s);
        double xi31=dum*(x1s*x2t - x1t*x2s);
        double xi32=dum*(x2r*x1t - x1r*x2t);
        double xi33=dum*(x1r*x2s - x2r*x1s);

        // get operator b of global derivatives
        for (int i=0; i<iel; i++)
        {
          int node_start = i*3;

          double hr   = deriv(0,i);
          double hs   = deriv(1,i);
          double ht   = deriv(2,i);

          double h1 = xi11*hr + xi12*hs + xi13*ht;
          double h2 = xi21*hr + xi22*hs + xi23*ht;
          double h3 = xi31*hr + xi32*hs + xi33*ht;

          bop(0,node_start+0) = h1 ;
          bop(0,node_start+1) = 0.0;
          bop(0,node_start+2) = 0.0;
          bop(1,node_start+0) = 0.0;
          bop(1,node_start+1) = h2 ;
          bop(1,node_start+2) = 0.0;
          bop(2,node_start+0) = 0.0;
          bop(2,node_start+1) = 0.0;
          bop(2,node_start+2) = h3 ;
          bop(3,node_start+0) = h2 ;
          bop(3,node_start+1) = h1 ;
          bop(3,node_start+2) = 0.0;
          bop(4,node_start+0) = 0.0;
          bop(4,node_start+1) = h3 ;
          bop(4,node_start+2) = h2 ;
          bop(5,node_start+0) = h3 ;
          bop(5,node_start+1) = 0.0;
          bop(5,node_start+2) = h1 ;
        }

        // call material law

        double ym  = material->m.stvenant->youngs;
        double pv  = material->m.stvenant->possionratio;

        // evaluate basic material values
        double d1=ym*(1.0 - pv)/((1.0 + pv)*(1.0 - 2.0*pv));
        double d2=ym*pv/((1.0 + pv)*(1.0 - 2.0*pv));
        double d3=ym/((1.0 + pv)*2.0);

        // set values in material-matrix
        D(0,0)=d1;
        D(0,1)=d2;
        D(0,2)=d2;
        D(0,3)=0.0;
        D(0,4)=0.0;
        D(0,5)=0.0;

        D(1,0)=d2;
        D(1,1)=d1;
        D(1,2)=d2;
        D(1,3)=0.0;
        D(1,4)=0.0;
        D(1,5)=0.0;

        D(2,0)=d2;
        D(2,1)=d2;
        D(2,2)=d1;
        D(2,3)=0.0;
        D(2,4)=0.0;
        D(2,5)=0.0;

        D(3,0)=0.0;
        D(3,1)=0.0;
        D(3,2)=0.0;
        D(3,3)=d3;
        D(3,4)=0.0;
        D(3,5)=0.0;

        D(4,0)=0.0;
        D(4,1)=0.0;
        D(4,2)=0.0;
        D(4,3)=0.0;
        D(4,4)=d3;
        D(4,5)=0.0;

        D(5,0)=0.0;
        D(5,1)=0.0;
        D(5,2)=0.0;
        D(5,3)=0.0;
        D(5,4)=0.0;
        D(5,5)=d3;

        // elastic stiffness matrix ke
        //ale3_keku(estif,bop,D,fac,nd);

        for (int j=0; j<nd; j++)
        {
          double db[6];
          for (int k=0; k<6; k++)
          {
            db[k] = 0.0;
            for (int l=0; l<6; l++)
            {
              db[k] += D(k,l)*bop(l,j)*fac ;
            }
          }
          for (int i=0; i<nd; i++)
          {
            double dum = 0.0;
            for (int m=0; m<6; m++)
            {
              dum = dum + bop(m,i)*db[m] ;
            }
            (*sys_mat)(i,j) += dum ;
          }
        }

        // hourglass stabalization stiffness matrix ke
        /*
          see also:
          (1) T. Belytschko and L.P. Bindeman:
              Assumed strain stabilization of the 8-node hexahedral element
              Comp. Meth. Appl. Mech. Eng.: 105 (1993) p. 225-260.
          (2) D.P. Flanagan and T. Belytschko:
              A uniform strain hexahedron and quadrilateral with orthogonal
              hourglass control
              Int. J. Num. Meth. Ing.: Vol. 17 (1981) p. 679-706.
         */
        if (Shape()==hex8 && nir == 1 && nis == 1 && nit == 1)
        {
          double ee = material->m.stvenant->youngs;
          double nu = material->m.stvenant->possionratio;
          double mu = ee / (2*(1+nu));

          // ASQBI
          double c1 = 1.0/(1.0 - nu);
          double c2 = (1.0 + nu)/3;
          double c3 = 1.0/(1.0 - nu);

          double         xc[3][8];
          double         b[3][8];

          double         a[3][3];
          double         ba[3];
          double         r[3][3];

          double         h[4][8] = {{1,1,-1,-1,-1,-1,1,1},{1,-1,-1,1,-1,1,1,-1},
                                    {1,-1,1,-1,1,-1,1,-1},{-1,1,-1,1,1,-1,1,-1}};
          double         lam[3][8] = {{-1,1,1,-1,-1,1,1,-1},
                                      {-1,-1,1,1,-1,-1,1,1},
                                      {-1,-1,-1,-1,1,1,1,1}};
          double         gam[4][8];
          double         lx[3];
          double         hh[3][3];

          double         gg00[8][8], gg11[8][8], gg22[8][8], gg33[8][8];
          double         gg01[8][8], gg10[8][8], gg02[8][8], gg20[8][8];
          double         gg12[8][8], gg21[8][8];

          double         kstab[24][24];

          // corotational coordinate system: rotation tensor r[3][3]
          // accord. to (1) Appendix A, Table 9
          for (int i=0; i<2; i++)
          {
            for (int j=0; j<3; j++)
            {
              a[i][j] = 0.0;
              for (int k=0; k<8; k++)
              {
                a[i][j] += lam[i][k]*xyze(j,k);
              }
            }
          }

          double dum =(a[0][0]*a[1][0]+a[0][1]*a[1][1]+a[0][2]*a[1][2])/
                      (a[0][0]*a[0][0]+a[0][1]*a[0][1]+a[0][2]*a[0][2]);

          a[1][0] = a[1][0] - dum * a[0][0];
          a[1][1] = a[1][1] - dum * a[0][1];
          a[1][2] = a[1][2] - dum * a[0][2];

          a[2][0] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
          a[2][1] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
          a[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

          for(int i = 0; i<3; i++)
          {
            ba[i] = sqrt(a[i][0]*a[i][0]+a[i][1]*a[i][1]+a[i][2]*a[i][2]);
            r[i][0] = a[i][0] / ba[i];
            r[i][1] = a[i][1] / ba[i];
            r[i][2] = a[i][2] / ba[i];
          }

          // transforming nodal coordinates to the corotational system
          for (int i=0; i<8; i++)
          {
            for (int j=0; j<3; j++)
            {
              xc[j][i] = r[j][0]*xyze(0,i)+r[j][1]*xyze(1,i)+r[j][2]*xyze(2,i);
            }
          }

          int l0=0,l1=0,l2=0,l3=0,l4=0,l5=0,l6=0,l7=0;

          // B-matrix b[3][8] accord. to (2) Appendix I, eqn (79)
          for (int i=0; i<3; i++)
          {
            int j = (i+1)%3;
            int k = (j+1)%3;
            for (int l=0; l<8;l++)
            {
              switch (l)
              {
              case 0:
                l0=0;l1=1;l2=2;l3=3;l4=4;l5=5;l6=6;l7=7;
                break;
              case 1:
                l0=1;l1=2;l2=3;l3=0;l4=5;l5=6;l6=7;l7=4;
                break;
              case 2:
                l0=2;l1=3;l2=0;l3=1;l4=6;l5=7;l6=4;l7=5;
                break;
              case 3:
                l0=3;l1=0;l2=1;l3=2;l4=7;l5=4;l6=5;l7=6;
                break;
              case 4:
                l0=4;l1=7;l2=6;l3=5;l4=0;l5=3;l6=2;l7=1;
                break;
              case 5:
                l0=5;l1=4;l2=7;l3=6;l4=1;l5=0;l6=3;l7=2;
                break;
              case 6:
                l0=6;l1=5;l2=4;l3=7;l4=2;l5=1;l6=0;l7=3;
                break;
              case 7:
                l0=7;l1=6;l2=5;l3=4;l4=3;l5=2;l6=1;l7=0;
                break;
              }
              b[i][l0] =1/(12 * vol) * (xc[j][l1]*((xc[k][l5] - xc[k][l2]) -
                                                   (xc[k][l3] - xc[k][l4])) +
                                        xc[j][l2]* (xc[k][l1] - xc[k][l3]) +
                                        xc[j][l3]*((xc[k][l2] - xc[k][l7]) -
                                                   (xc[k][l4] - xc[k][l1])) +
                                        xc[j][l4]*((xc[k][l7] - xc[k][l5]) -
                                                   (xc[k][l1] - xc[k][l3])) +
                                        xc[j][l5]* (xc[k][l4] - xc[k][l1]) +
                                        xc[j][l7]* (xc[k][l3] - xc[k][l4]) );
            }
          }

          // gamma vectors, accord. to (1) eqn (2.12b)
          for (int i=0; i<4; i++)
          {
            for (int j=0; j<8; j++)
            {
              gam[i][j] = 0.125 * h[i][j];
              for (int k=0; k<3; k++)
              {
                dum = h[i][0]*xc[k][0]+h[i][1]*xc[k][1]+h[i][2]*xc[k][2]+
                      h[i][3]*xc[k][3]+h[i][4]*xc[k][4]+h[i][5]*xc[k][5]+
                      h[i][6]*xc[k][6]+h[i][7]*xc[k][7];
                gam[i][j] -= 0.125 * dum * b[k][j];
              }
            }
          }

          /* lambda * x (auxiliary vector) */
          lx[0] = lam[0][0]*xc[0][0]+lam[0][1]*xc[0][1]+lam[0][2]*xc[0][2]+
                  lam[0][3]*xc[0][3]+lam[0][4]*xc[0][4]+lam[0][5]*xc[0][5]+
                  lam[0][6]*xc[0][6]+lam[0][7]*xc[0][7];
          lx[1] = lam[1][0]*xc[1][0]+lam[1][1]*xc[1][1]+lam[1][2]*xc[1][2]+
                  lam[1][3]*xc[1][3]+lam[1][4]*xc[1][4]+lam[1][5]*xc[1][5]+
                  lam[1][6]*xc[1][6]+lam[1][7]*xc[1][7];
          lx[2] = lam[2][0]*xc[2][0]+lam[2][1]*xc[2][1]+lam[2][2]*xc[2][2]+
                  lam[2][3]*xc[2][3]+lam[2][4]*xc[2][4]+lam[2][5]*xc[2][5]+
                  lam[2][6]*xc[2][6]+lam[2][7]*xc[2][7];

          /* H_ij, accord. to (1) eqns. (3.15d) and (3.15e) */
          hh[0][0] = 1.0/3.0 * (lx[1]*lx[2])/lx[0];
          hh[1][1] = 1.0/3.0 * (lx[2]*lx[0])/lx[1];
          hh[2][2] = 1.0/3.0 * (lx[0]*lx[1])/lx[2];
          hh[0][1] = 1.0/3.0 * lx[2];
          hh[1][0] = 1.0/3.0 * lx[2];
          hh[0][2] = 1.0/3.0 * lx[1];
          hh[2][0] = 1.0/3.0 * lx[1];
          hh[1][2] = 1.0/3.0 * lx[0];
          hh[2][1] = 1.0/3.0 * lx[0];

          /* stabalization matrix with respect to the corotational ccord. system. */
          /* rearranging orders of dofs, accord. to (1) eqns. (3.15a) to (3.15c) */
          for (int i=0; i<8; i++)
          {
            for (int j=0; j<8; j++)
            {
              gg00[i][j] = gam[0][i] * gam[0][j];
              gg11[i][j] = gam[1][i] * gam[1][j];
              gg22[i][j] = gam[2][i] * gam[2][j];
              gg33[i][j] = gam[3][i] * gam[3][j];
              gg01[i][j] = gam[0][i] * gam[1][j];
              gg10[i][j] = gam[1][i] * gam[0][j];
              gg02[i][j] = gam[0][i] * gam[2][j];
              gg20[i][j] = gam[2][i] * gam[0][j];
              gg12[i][j] = gam[1][i] * gam[2][j];
              gg21[i][j] = gam[2][i] * gam[1][j];

              /* kstab 00 */
              kstab[i*3][j*3]     = 2*mu* (hh[0][0]*(c1*(gg11[i][j] + gg22[i][j])
                                                     + c2*gg33[i][j]) + 0.5 * (hh[1][1] + hh[2][2]) * gg00[i][j]);

              /* kstab 11 */
              kstab[i*3+1][j*3+1] = 2*mu* (hh[1][1]*(c1*(gg22[i][j] + gg00[i][j])
                                                     + c2*gg33[i][j]) + 0.5 * (hh[2][2] + hh[0][0]) * gg11[i][j]);

              /* kstab 22 */
              kstab[i*3+2][j*3+2] = 2*mu* (hh[2][2]*(c1*(gg00[i][j] + gg11[i][j])
                                                     + c2*gg33[i][j]) + 0.5 * (hh[0][0] + hh[1][1]) * gg22[i][j]);

              /* kstab 01 */
              kstab[i*3][j*3+1]   = 2*mu* (hh[0][1]*(c3*gg10[i][j]+0.5*gg01[i][j]));

              /* kstab 10 */
              kstab[i*3+1][j*3]   = 2*mu* (hh[1][0]*(c3*gg01[i][j]+0.5*gg10[i][j]));

              /* kstab 02 */
              kstab[i*3][j*3+2]   = 2*mu* (hh[0][2]*(c3*gg20[i][j]+0.5*gg02[i][j]));

              /* kstab 20 */
              kstab[i*3+2][j*3]   = 2*mu* (hh[2][0]*(c3*gg02[i][j]+0.5*gg20[i][j]));

              /* kstab 12 */
              kstab[i*3+1][j*3+2] = 2*mu* (hh[1][2]*(c3*gg21[i][j]+0.5*gg12[i][j]));

              /* kstab 21 */
              kstab[i*3+2][j*3+1] = 2*mu* (hh[2][1]*(c3*gg12[i][j]+0.5*gg21[i][j]));
            }
          }

          /* transforming kstab to the global coordinate system and */
          /* and adding to the one point quadrature matrix */
          for (int i=0; i<8; i++)
          {
            for (int j=0;j<8;j++)
            {
              for (int k=0;k<3;k++)
              {
                for (int l=0;l<3;l++)
                {
                  (*sys_mat)(i*3+k,j*3+l) += (r[0][k]*kstab[i*3+0][j*3+0] +
                                              r[1][k]*kstab[i*3+1][j*3+0] +
                                              r[2][k]*kstab[i*3+2][j*3+0]) * r[0][l] +
                                             (r[0][k]*kstab[i*3+0][j*3+1] +
                                              r[1][k]*kstab[i*3+1][j*3+1] +
                                              r[2][k]*kstab[i*3+2][j*3+1]) * r[1][l] +
                                             (r[0][k]*kstab[i*3+0][j*3+2] +
                                              r[1][k]*kstab[i*3+1][j*3+2] +
                                              r[2][k]*kstab[i*3+2][j*3+2]) * r[2][l];
                }
              }
            }
          }
        }
      }
    }
  }

}



void DRT::Elements::Ale3::intg(ALE3_DATA& data)
{
  DOUBLE  q14, q16, q124;
  DOUBLE  palpha,pbeta;

  switch (Shape())
  {
  case hex8:
  case hex20:
  case hex27:
    /*----------------------------------------------------------------------*
     |     INTEGRATION PARAMETERS FOR    H E X A H E D R A L     ELEMENTS   |
     |     GAUSS SAMPLING POINTS  AT     R/S/T-COORDINATES   RESPECTIVELY   |
     |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
     *----------------------------------------------------------------------*/
    switch (ngp_[0])
    {
    case 3:
      data.xgpr[0] = -0.7745966692415;
      data.xgpr[1] =  0.0;
      data.xgpr[2] =  0.7745966692415;

      data.wgtr[0] =  0.5555555555556;
      data.wgtr[1] =  0.8888888888889;
      data.wgtr[2] =  0.5555555555556;
      break;
    case 2:
      data.xgpr[0] = -0.5773502691896;
      data.xgpr[1] =  0.5773502691896;

      data.wgtr[0] = 1.0            ;
      data.wgtr[1] = 1.0            ;
      break;
    case 1:
      data.xgpr[0] = 0.0;

      data.wgtr[0] = 2.0;
      break;
    default:
      dserror("unknown number of gaussian points in ale3_intg");
      break;
    }
    switch (ngp_[1])
    {
    case 3:
      data.xgps[0] = -0.7745966692415;
      data.xgps[1] =  0.0;
      data.xgps[2] =  0.7745966692415;

      data.wgts[0] =  0.5555555555556;
      data.wgts[1] =  0.8888888888889;
      data.wgts[2] =  0.5555555555556;
      break;
    case 2:
      data.xgps[0] = -0.5773502691896;
      data.xgps[1] =  0.5773502691896;

      data.wgts[0] = 1.0            ;
      data.wgts[1] = 1.0            ;
      break;
    case 1:
      data.xgps[0] = 0.0;

      data.wgts[0] = 2.0;
      break;
    default:
      dserror("unknown number of gaussian points in ale3_intg");
      break;
    }
    switch (ngp_[2])
    {
    case 3:
      data.xgpt[0] = -0.7745966692415;
      data.xgpt[1] =  0.0;
      data.xgpt[2] =  0.7745966692415;

      data.wgtt[0] =  0.5555555555556;
      data.wgtt[1] =  0.8888888888889;
      data.wgtt[2] =  0.5555555555556;
      break;
    case 2:
      data.xgpt[0] = -0.5773502691896;
      data.xgpt[1] =  0.5773502691896;

      data.wgtt[0] = 1.0            ;
      data.wgtt[1] = 1.0            ;
      break;
    case 1:
      data.xgpt[0] = 0.0;

      data.wgtt[0] = 2.0;
      break;
    default:
      dserror("unknown number of gaussian points in ale3_intg");
      break;
    }
    break;
  case tet4:
  case tet10:
    q14 = 1.0/4.0;
    q16 = 1.0/6.0;
    q124= 1.0/24.0;
    palpha = (5.0+3.0*sqrt(5.0))/20.0;
    pbeta  = (5.0-sqrt(5.0))/20.0;

    /*----------------------------------------------------------------------*
     |     INTEGRATION PARAMETERS FOR    T E T R A H E D R A L   ELEMENTS   |
     |     GAUSS SAMPLING POINTS  AT     R/S/T-COORDINATES   RESPECTIVELY   |
     |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
     *----------------------------------------------------------------------*/

    /*----------------------------------------------------------------------*
     |    GAUSS INTEGRATION         1 SAMPLING POINT, DEG.OF PRECISION 1    |
     |                              CASE 0                                  |
     *----------------------------------------------------------------------*/
    switch (ngp_[0])
    {
    case 1:
      data.xgpr[0]    =  q14 ;
      data.xgps[0]    =  q14 ;
      data.xgpt[0]    =  q14 ;
      data.wgtr[0]    =  q16 ;
      break;
      /*----------------------------------------------------------------*
       | GAUSS INTEGRATION    4 SAMPLING POINTS, DEG.OF PRECISION 2    |
       |                      CASE 1                                   |
       *----------------------------------------------------------------*/
    case 4:
      data.xgpr[0]    =    pbeta ;
      data.xgpr[1]    =    palpha;
      data.xgpr[2]    =    pbeta ;
      data.xgpr[3]    =    pbeta ;
      data.xgps[0]    =    pbeta ;
      data.xgps[1]    =    pbeta ;
      data.xgps[2]    =    palpha;
      data.xgps[3]    =    pbeta ;
      data.xgpt[0]    =    pbeta ;
      data.xgpt[1]    =    pbeta ;
      data.xgpt[2]    =    pbeta ;
      data.xgpt[3]    =    palpha;
      data.wgtr[0]    =    q124  ;
      data.wgtr[1]    =    q124  ;
      data.wgtr[2]    =    q124  ;
      data.wgtr[3]    =    q124  ;
      break;
    default:
      dserror("unknown number of gausian points");
      break;
    }
    break;
  default:
    dserror("unknown typ of discretisation");
    break;
  }
}


void DRT::Elements::Ale3::funct_deriv(
    vector<double>& funct,
    Epetra_SerialDenseMatrix& deriv,
    double      r,
    double      s,
    double      t,
    DiscretizationType     typ,
    int         option
    )
{
  const DOUBLE   q18 = 1.0/8.0;
  DOUBLE         rp,sp,tp,rm,sm,tm,rrm,ssm,ttm;
  DOUBLE         t1,t2,t3,t4;

  /* if option ==0 only funtion evaluation, if option==1 also derivatives */
  rp  = 1.0+r;
  rm  = 1.0-r;
  sp  = 1.0+s;
  sm  = 1.0-s;
  tp  = 1.0+t;
  tm  = 1.0-t;
  rrm = 1.0-r*r;
  ssm = 1.0-s*s;
  ttm = 1.0-t*t;

  switch (typ)
  {
    /*  LINEAR SHAPE FUNCTIONS  SERENDIPITY (8-NODED ELEMENT) */
    case hex8:
      funct[0]=q18*rp*sm*tm;
      funct[1]=q18*rp*sp*tm;
      funct[2]=q18*rm*sp*tm;
      funct[3]=q18*rm*sm*tm;
      funct[4]=q18*rp*sm*tp;
      funct[5]=q18*rp*sp*tp;
      funct[6]=q18*rm*sp*tp;
      funct[7]=q18*rm*sm*tp;

      if (option==1)  /* check for derivative evaluation */
      {
        deriv(0,0)= q18*sm*tm  ;
        deriv(0,1)= q18*sp*tm  ;
        deriv(0,2)=-deriv(0,1);
        deriv(0,3)=-deriv(0,0);
        deriv(0,4)= q18*sm*tp  ;
        deriv(0,5)= q18*sp*tp  ;
        deriv(0,6)=-deriv(0,5);
        deriv(0,7)=-deriv(0,4);
        deriv(1,0)=-q18*tm*rp  ;
        deriv(1,1)=-deriv(1,0);
        deriv(1,2)= q18*tm*rm  ;
        deriv(1,3)=-deriv(1,2);
        deriv(1,4)=-q18*tp*rp  ;
        deriv(1,5)=-deriv(1,4);
        deriv(1,6)= q18*tp*rm  ;
        deriv(1,7)=-deriv(1,6);
        deriv(2,0)=-q18*rp*sm  ;
        deriv(2,1)=-q18*rp*sp  ;
        deriv(2,2)=-q18*rm*sp  ;
        deriv(2,3)=-q18*rm*sm  ;
        deriv(2,4)=-deriv(2,0);
        deriv(2,5)=-deriv(2,1);
        deriv(2,6)=-deriv(2,2);
        deriv(2,7)=-deriv(2,3);
      }
      break;


    case hex20: /* QUADRATIC shape functions and their natural derivatives
                   without central nodes */


      /*   Shape functions and their derivatives for a 20 noded hexaedron
       *   ==============================================================
       *
       *   Numbering of the nodes:
       *   -----------------------
       *   - this is the numbering used in GiD!!
       *   - the numbering of the brick1 element is different!!
       *
       *
       *
       *                          ^ t          / s
       *                          |           /
       *                          |          /
       *                    8     |   19    /   7
       *                    o-----|---o---------o
       *                   /|     |       /    /|
       *                  / |     |      /    / |
       *                 /  |     |     /    /  |
       *              20o   |     |    /    o18 |
       *               /  16o     |   /    /    o15
       *              /     |     |  /    /     |
       *             /      |  17 | /  6 /      |
       *          5 o---------o---------o       |
       *            |       |     *-----|---------------->
       *            |       o---------o-|-------o         r
       *            |      / 4       11 |      /3
       *            |     /             |     /
       *          13o    /              o14  /
       *            | 12o               |   o10
       *            |  /                |  /
       *            | /                 | /
       *            |/                  |/
       *            o---------o---------o
       *           1         9         2
       */



      /* form basic values */
      rp  = 1.+r;
      rm  = 1.-r;
      sp  = 1.+s;
      sm  = 1.-s;
      tp  = 1.+t;
      tm  = 1.-t;
      rrm = 1.-r*r;
      ssm = 1.-s*s;
      ttm = 1.-t*t;

      funct[ 0] = -1./8.*rm*sm*tm*(2.+r+s+t);
      funct[ 1] = -1./8.*rp*sm*tm*(2.-r+s+t);
      funct[ 2] = -1./8.*rp*sp*tm*(2.-r-s+t);
      funct[ 3] = -1./8.*rm*sp*tm*(2.+r-s+t);
      funct[ 4] = -1./8.*rm*sm*tp*(2.+r+s-t);
      funct[ 5] = -1./8.*rp*sm*tp*(2.-r+s-t);
      funct[ 6] = -1./8.*rp*sp*tp*(2.-r-s-t);
      funct[ 7] = -1./8.*rm*sp*tp*(2.+r-s-t);

      funct[ 8] =  1./4.*rrm*sm*tm;
      funct[ 9] =  1./4.*rp*ssm*tm;
      funct[10] =  1./4.*rrm*sp*tm;
      funct[11] =  1./4.*rm*ssm*tm;

      funct[12] =  1./4.*rm*sm*ttm;
      funct[13] =  1./4.*rp*sm*ttm;
      funct[14] =  1./4.*rp*sp*ttm;
      funct[15] =  1./4.*rm*sp*ttm;

      funct[16] =  1./4.*rrm*sm*tp;
      funct[17] =  1./4.*rp*ssm*tp;
      funct[18] =  1./4.*rrm*sp*tp;
      funct[19] =  1./4.*rm*ssm*tp;


      /* first derivative evaluation */
      deriv(0, 0) =  1./8.*   sm*tm*(2.*r+s+t+1.);
      deriv(1, 0) =  1./8.*rm*   tm*(r+2.*s+t+1.);
      deriv(2, 0) =  1./8.*rm*sm*   (r+s+2.*t+1.);

      deriv(0, 1) =  1./8.*   sm*tm*(2.*r-s-t-1.);
      deriv(1, 1) = -1./8.*rp*   tm*(r-2.*s-t-1.);
      deriv(2, 1) = -1./8.*rp*sm*   (r-s-2.*t-1.);

      deriv(0, 2) =  1./8.*   sp*tm*(2.*r+s-t-1.);
      deriv(1, 2) =  1./8.*rp*   tm*(r+2.*s-t-1.);
      deriv(2, 2) = -1./8.*rp*sp*   (r+s-2.*t-1.);

      deriv(0, 3) =  1./8.*   sp*tm*(2.*r-s+t+1.);
      deriv(1, 3) = -1./8.*rm*   tm*(r-2.*s+t+1.);
      deriv(2, 3) =  1./8.*rm*sp*   (r-s+2.*t+1.);

      deriv(0, 4) =  1./8.*   sm*tp*(2.*r+s-t+1.);
      deriv(1, 4) =  1./8.*rm*   tp*(r+2.*s-t+1.);
      deriv(2, 4) = -1./8.*rm*sm*   (r+s-2.*t+1.);

      deriv(0, 5) =  1./8.*   sm*tp*(2.*r-s+t-1.);
      deriv(1, 5) = -1./8.*rp*   tp*(r-2.*s+t-1.);
      deriv(2, 5) =  1./8.*rp*sm*   (r-s+2.*t-1.);

      deriv(0, 6) =  1./8.*   sp*tp*(2.*r+s+t-1.);
      deriv(1, 6) =  1./8.*rp*   tp*(r+2.*s+t-1.);
      deriv(2, 6) =  1./8.*rp*sp*   (r+s+2.*t-1.);

      deriv(0, 7) =  1./8.*   sp*tp*(2.*r-s-t+1.);
      deriv(1, 7) = -1./8.*rm*   tp*(r-2.*s-t+1.);
      deriv(2, 7) = -1./8.*rm*sp*   (r-s-2.*t+1.);


      deriv(0, 8) = -1./2.*r*sm*tm;
      deriv(1, 8) = -1./4.*rm*rp*tm;
      deriv(2, 8) = -1./4.*rm*rp*sm;

      deriv(0, 9) =  1./4.*sm*sp*tm;
      deriv(1, 9) = -1./2.*rp*s*tm;
      deriv(2, 9) = -1./4.*rp*sm*sp;

      deriv(0,10) = -1./2.*r*sp*tm;
      deriv(1,10) =  1./4.*rm*rp*tm;
      deriv(2,10) = -1./4.*rm*rp*sp;

      deriv(0,11) = -1./4.*sm*sp*tm;
      deriv(1,11) = -1./2.*s*tm*rm;
      deriv(2,11) = -1./4.*sm*sp*rm;


      deriv(0,12) = -1./4.*sm*tm*tp;
      deriv(1,12) = -1./4.*rm*tm*tp;
      deriv(2,12) = -1./2.*t*rm*sm;

      deriv(0,13) =  1./4.*sm*tm*tp;
      deriv(1,13) = -1./4.*rp*tm*tp;
      deriv(2,13) = -1./2.*t*rp*sm;

      deriv(0,14) =  1./4.*sp*tm*tp;
      deriv(1,14) =  1./4.*rp*tm*tp;
      deriv(2,14) = -1./2.*t*rp*sp;

      deriv(0,15) = -1./4.*sp*tm*tp;
      deriv(1,15) =  1./4.*rm*tm*tp;
      deriv(2,15) = -1./2.*t*rm*sp;


      deriv(0,16) = -1./2.*r*sm*tp;
      deriv(1,16) = -1./4.*rm*rp*tp;
      deriv(2,16) =  1./4.*rm*rp*sm;

      deriv(0,17) =  1./4.*sm*sp*tp;
      deriv(1,17) = -1./2.*s*tp*rp;
      deriv(2,17) =  1./4.*sm*sp*rp;

      deriv(0,18) = -1./2.*r*sp*tp;
      deriv(1,18) =  1./4.*rm*rp*tp;
      deriv(2,18) =  1./4.*rm*rp*sp;

      deriv(0,19) = -1./4.*sm*sp*tp;
      deriv(1,19) = -1./2.*s*tp*rm;
      deriv(2,19) =  1./4.*sm*sp*rm;
      break;

  case hex27: /* QUADRATIC shape functions and their natural derivatives
                 with central nodes                         ----*/
  {
    /* form basic values */
    DOUBLE drm1,dr00,drp1,dsm1,ds00,dsp1,dtm1,dt00,dtp1;
    DOUBLE rm1,r00,rp1,sm1,s00,sp1,tm1,t00,tp1;

    rm1=1./2.*r*(r - 1.);
    r00=(1. - r*r);
    rp1=1./2.*r*(r + 1.);
    sm1=1./2.*s*(s - 1.);
    s00=(1. - s*s);
    sp1=1./2.*s*(s + 1.);
    tm1=1./2.*t*(t - 1.);
    t00=(1. - t*t);
    tp1=1./2.*t*(t + 1.);

    drm1 = r - 1./2.;
    dr00 = -2. * r;
    drp1 = r + 1./2.;
    dsm1 = s - 1./2.;
    ds00 = -2. * s;
    dsp1 = s + 1./2.;
    dtm1 = t - 1./2.;
    dt00 = -2. * t;
    dtp1 = t + 1./2.;

    funct[0] = rp1*sp1*tp1;
    funct[1] = sm1*rp1*tp1;
    funct[2] = rm1*sm1*tp1;
    funct[3] = rm1*sp1*tp1;
    funct[4] = tm1*rp1*sp1;
    funct[5] = sm1*tm1*rp1;
    funct[6] = rm1*sm1*tm1;
    funct[7] = rm1*tm1*sp1;
    funct[8] = s00*rp1*tp1;
    funct[9] = r00*sm1*tp1;
    funct[10] = s00*rm1*tp1;
    funct[11] = r00*sp1*tp1;
    funct[12] = t00*rp1*sp1;
    funct[13] = t00*sm1*rp1;
    funct[14] = t00*rm1*sm1;
    funct[15] = t00*rm1*sp1;
    funct[16] = s00*tm1*rp1;
    funct[17] = r00*sm1*tm1;
    funct[18] = s00*rm1*tm1;
    funct[19] = r00*tm1*sp1;
    funct[20] = r00*s00*tp1;
    funct[21] = s00*t00*rp1;
    funct[22] = r00*t00*sm1;
    funct[23] = s00*t00*rm1;
    funct[24] = r00*t00*sp1;
    funct[25] = r00*s00*tm1;
    funct[26] = r00*s00*t00;

    if (option==1) /* --> first derivative evaluation */
    {
      deriv(0,0) = sp1*tp1*drp1;
      deriv(0,1) = sm1*tp1*drp1;
      deriv(0,2) = sm1*tp1*drm1;
      deriv(0,3) = sp1*tp1*drm1;
      deriv(0,4) = tm1*sp1*drp1;
      deriv(0,5) = sm1*tm1*drp1;
      deriv(0,6) = sm1*tm1*drm1;
      deriv(0,7) = tm1*sp1*drm1;
      deriv(0,8) = s00*tp1*drp1;
      deriv(0,9) = sm1*tp1*dr00;
      deriv(0,10) = s00*tp1*drm1;
      deriv(0,11) = sp1*tp1*dr00;
      deriv(0,12) = t00*sp1*drp1;
      deriv(0,13) = t00*sm1*drp1;
      deriv(0,14) = t00*sm1*drm1;
      deriv(0,15) = t00*sp1*drm1;
      deriv(0,16) = s00*tm1*drp1;
      deriv(0,17) = sm1*tm1*dr00;
      deriv(0,18) = s00*tm1*drm1;
      deriv(0,19) = tm1*sp1*dr00;
      deriv(0,20) = s00*tp1*dr00;
      deriv(0,21) = s00*t00*drp1;
      deriv(0,22) = t00*sm1*dr00;
      deriv(0,23) = s00*t00*drm1;
      deriv(0,24) = t00*sp1*dr00;
      deriv(0,25) = s00*tm1*dr00;
      deriv(0,26) = s00*t00*dr00;

      deriv(1,0) = rp1*tp1*dsp1;
      deriv(1,1) = rp1*tp1*dsm1;
      deriv(1,2) = rm1*tp1*dsm1;
      deriv(1,3) = rm1*tp1*dsp1;
      deriv(1,4) = tm1*rp1*dsp1;
      deriv(1,5) = tm1*rp1*dsm1;
      deriv(1,6) = rm1*tm1*dsm1;
      deriv(1,7) = rm1*tm1*dsp1;
      deriv(1,8) = rp1*tp1*ds00;
      deriv(1,9) = r00*tp1*dsm1;
      deriv(1,10) = rm1*tp1*ds00;
      deriv(1,11) = r00*tp1*dsp1;
      deriv(1,12) = t00*rp1*dsp1;
      deriv(1,13) = t00*rp1*dsm1;
      deriv(1,14) = t00*rm1*dsm1;
      deriv(1,15) = t00*rm1*dsp1;
      deriv(1,16) = tm1*rp1*ds00;
      deriv(1,17) = r00*tm1*dsm1;
      deriv(1,18) = rm1*tm1*ds00;
      deriv(1,19) = r00*tm1*dsp1;
      deriv(1,20) = r00*tp1*ds00;
      deriv(1,21) = t00*rp1*ds00;
      deriv(1,22) = r00*t00*dsm1;
      deriv(1,23) = t00*rm1*ds00;
      deriv(1,24) = r00*t00*dsp1;
      deriv(1,25) = r00*tm1*ds00;
      deriv(1,26) = r00*t00*ds00;

      deriv(2,0) = rp1*sp1*dtp1;
      deriv(2,1) = sm1*rp1*dtp1;
      deriv(2,2) = rm1*sm1*dtp1;
      deriv(2,3) = rm1*sp1*dtp1;
      deriv(2,4) = rp1*sp1*dtm1;
      deriv(2,5) = sm1*rp1*dtm1;
      deriv(2,6) = rm1*sm1*dtm1;
      deriv(2,7) = rm1*sp1*dtm1;
      deriv(2,8) = s00*rp1*dtp1;
      deriv(2,9) = r00*sm1*dtp1;
      deriv(2,10) = s00*rm1*dtp1;
      deriv(2,11) = r00*sp1*dtp1;
      deriv(2,12) = rp1*sp1*dt00;
      deriv(2,13) = sm1*rp1*dt00;
      deriv(2,14) = rm1*sm1*dt00;
      deriv(2,15) = rm1*sp1*dt00;
      deriv(2,16) = s00*rp1*dtm1;
      deriv(2,17) = r00*sm1*dtm1;
      deriv(2,18) = s00*rm1*dtm1;
      deriv(2,19) = r00*sp1*dtm1;
      deriv(2,20) = r00*s00*dtp1;
      deriv(2,21) = s00*rp1*dt00;
      deriv(2,22) = r00*sm1*dt00;
      deriv(2,23) = s00*rm1*dt00;
      deriv(2,24) = r00*sp1*dt00;
      deriv(2,25) = r00*s00*dtm1;
      deriv(2,26) = r00*s00*dt00;
    }
    break;
  }

    case tet4: /* LINEAR SHAPE FUNCTIONS */
/*
      t1=r;
      t2=s;
      t3=t;
      t4=1.-r-s-t;
*/
   t1=1.-r-s-t;
   t2=r;
   t3=s;
   t4=t;

      funct[0]= t1;
      funct[1]= t2;
      funct[2]= t3;
      funct[3]= t4;

      if(option==1) /* first derivative evaluation */
      {
        /*
        deriv(0,0)= 1.;
        deriv(0,1)= ZERO;
        deriv(0,2)= ZERO;
        deriv(0,3)=-1.;

        deriv(1,0)= ZERO;
        deriv(1,1)= 1.;
        deriv(1,2)= ZERO;
        deriv(1,3)=-1.;

        deriv(2,0)= ZERO;
        deriv(2,1)= ZERO;
        deriv(2,2)= 1.;
        deriv(2,3)=-1.;
        */
      deriv(0,0)=-1.;
      deriv(0,1)= 1.;
      deriv(0,2)= ZERO;
      deriv(0,3)= ZERO;

      deriv(1,0)=-1.;
      deriv(1,1)= ZERO;
      deriv(1,2)= 1.;
      deriv(1,3)= ZERO;

      deriv(2,0)=-1.;
      deriv(2,1)= ZERO;
      deriv(2,2)= ZERO;
      deriv(2,3)= 1.;

      } /* endif (option==1) */
      break;

    case tet10: /*  QUADRATIC SHAPE FUNCTIONS */

      dserror("shape functions for tet10 not yet implemented \n");
      /* form basic values */
#if 0
      t1=r;
      t2=s;
      t3=t;
      t4=1.-r-s-t;

      funct[0] =  ;
      funct[1] =  ;
      funct[2] = ;
      funct[3] = ;
      funct[4] = ;
      funct[5] = ;
      funct[6] = ;
      funct[7] = ;
      funct[8] = ;
      funct[9] = ;


      if(option==1) /* first derivative evaluation */
      {
        deriv[0][0] = ;
        deriv[1][0] = ;
        deriv[2][0] = ;

        deriv[0][1] = ;
        deriv[1][1] = ;
        deriv[2][1] = ;

        deriv[0][2] = ;
        deriv[1][2] = ;
        deriv[2][2] = ;

        deriv[0][3] = ;
        deriv[1][3] = ;
        deriv[2][3] = ;

        deriv[0][4] = ;
        deriv[1][4] = ;
        deriv[2][4] = ;

        deriv[0][5] = ;
        deriv[1][5] = ;
        deriv[2][5] = ;

        deriv[0][6] = ;
        deriv[1][6] = ;
        deriv[2][6] = ;

        deriv[0][7] = ;
        deriv[1][7] = ;
        deriv[2][7] = ;

        deriv[0][8] = ;
        deriv[1][8] = ;
        deriv[2][8] = ;

        deriv[0][9] = ;
        deriv[1][9] = ;
        deriv[2][9] = ;
      }
      break;
#endif

    default:
      dserror("unknown typ of interpolation");
      break;
  }
}


//=======================================================================
//=======================================================================

int DRT::Elements::Ale3Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif
#endif
#endif
