#ifdef D_ALE
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include <mpi.h>
#endif

#include "ale3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_mat/stvenantkirchhoff.H"

using namespace DRT::Utils;

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
  RefCountPtr<MAT::Material> mat = Material();

  switch(act)
  {
  case calc_ale_lin_stiff:
  {
    //RefCountPtr<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
    //vector<double> my_dispnp(lm.size());
    //DRT::Utils::ExtractMyValues(*dispnp,my_dispnp,lm);

    static_ke(lm,&elemat1,&elevec1,mat,params);

    break;
  }
  default:
    dserror("Unknown type of action for Ale3");
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
				    RefCountPtr<MAT::Material> material,
                                    ParameterList&            params)
{
  const int iel = NumNode();
  const int nd  = 3 * iel;
  const DiscretizationType distype = this->Shape();


  //  get material using class StVenantKirchhoff
  if (material->MaterialType()!=m_stvenant)
    dserror("stvenant material expected but got type %d", material->MaterialType());
  MAT::StVenantKirchhoff* actmat = static_cast<MAT::StVenantKirchhoff*>(material.get());

  Epetra_SerialDenseMatrix xyze(3,iel);

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
    xyze(2,i)=Nodes()[i]->X()[2];
  }

  /*----------------------------------------- declaration of variables ---*/
  Epetra_SerialDenseVector  funct(iel);
  Epetra_SerialDenseMatrix 	deriv(3,iel);
  Epetra_SerialDenseMatrix 	xjm(3,3);
  Epetra_SerialDenseMatrix 	bop(6,3*iel);
  Epetra_SerialDenseMatrix 	D(6,6);

  double                        vol=0.;

  // gaussian points
  const GaussRule3D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints3D  intpoints = getIntegrationPoints3D(gaussrule);

  // integration loops
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
        const double e1 = intpoints.qxg[iquad][0];
        const double e2 = intpoints.qxg[iquad][1];
        const double e3 = intpoints.qxg[iquad][2];
        // shape functions and their derivatives
        DRT::Utils::shape_function_3D(funct,e1,e2,e3,distype);
        DRT::Utils::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

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
        const double fac = intpoints.qwgt[iquad]*det;
        vol += fac;

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

        const double dum=1.0/det;

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
          const int node_start = i*3;

          const double hr   = deriv(0,i);
          const double hs   = deriv(1,i);
          const double ht   = deriv(2,i);

          const double h1 = xi11*hr + xi12*hs + xi13*ht;
          const double h2 = xi21*hr + xi22*hs + xi23*ht;
          const double h3 = xi31*hr + xi32*hs + xi33*ht;

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
	actmat->SetupCmat(&D);

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

	//Integration rule for hour-glass-stabilization. Not used in the moment. If needed, 
	//it should be implemented within the getOptimalGaussrule-method  
#if 0
        if (distype==hex8 && intpoints.nquad == 1)
        {
          const double ee = material->m.stvenant->youngs;
          const double nu = material->m.stvenant->possionratio;
          const double mu = ee / (2*(1+nu));

          // ASQBI
          const double c1 = 1.0/(1.0 - nu);
          const double c2 = (1.0 + nu)/3;
          const double c3 = 1.0/(1.0 - nu);

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

          const double dum =(a[0][0]*a[1][0]+a[0][1]*a[1][1]+a[0][2]*a[1][2])/
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
            const int j = (i+1)%3;
            const int k = (j+1)%3;
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
                const double dum = h[i][0]*xc[k][0]+h[i][1]*xc[k][1]+h[i][2]*xc[k][2]+
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
#endif
  }

}



// get optimal gaussrule for discretization type
GaussRule3D DRT::Elements::Ale3::getOptimalGaussrule(const DiscretizationType& distype)
{
    GaussRule3D rule = intrule3D_undefined;
    switch (distype)
    {
    case hex8:
        rule = intrule_hex_8point;
        break;
    case hex20: case hex27:
        rule = intrule_hex_27point;
        break;
    case tet4:
        rule = intrule_tet_4point;
        break;
    case tet10:
        rule = intrule_tet_5point;
        break;
    default: 
        dserror("unknown number of nodes for gaussrule initialization");
    }
    return rule;
}


//=======================================================================
//=======================================================================

int DRT::Elements::Ale3Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif
#endif
