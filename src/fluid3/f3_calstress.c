/*!----------------------------------------------------------------------
  \file
  \brief calculation of fluid3 stresses

  <pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
  </pre>

------------------------------------------------------------------------*/

/*! 
  \addtogroup FLUID3
  *//*! @{ (documentation module open)*/
#ifdef D_FLUID3
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"
/*----------------------------------------------------------------------*
  |                                                     m.gee 06/01    |
  | vector of material laws                                            |
  | defined in global_control.c                                        |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!--------------------------------------------------------------------- 
  \brief calculation of fluid stresses for fluid3

  <pre>                                                       mn 03/04

  calculation of the stress tensor sigma, which es stored at the nodes

  sigma = -p_real*I + 2*nue * eps(u) 				      

  </pre>

  \param    viscstr    INT         (i)    include viscose stresses yes/no
  \param   *data       FLUID_DATA  (i)
  \param   *ele        ELMENT      (i/o)  actual element
  \param  **evel       DOUBLE      (-)    element velocities
  \param   *epre       DOUBLE      (-)    element pressure
  \param  **funct      DOUBLE      (-)    shape functions
  \param  **deriv      DOUBLE      (-)    natural deriv. of shape funct.
  \param  **derxy      DOUBLE      (-)    global deriv. of sape funct.       
  \param  **vderxy     DOUBLE      (-)    global vel. deriv
  \param  **xjm        DOUBLE      (-)    jacobian matrix
  \param  **wa1        DOUBLE      (-)    working array
  \param  **xyze       DOUBLE      (-)    element coordinates
  \param  **sigmaint   DOUBLE      (-)    stresses at GAUSS point 

  \return void                                               

  ------------------------------------------------------------------------*/
void f3_calelestress(
    INT             viscstr,
    FLUID_DATA     *data, 
    ELEMENT        *ele,
    DOUBLE        **evel, 
    DOUBLE         *epre,
    DOUBLE         *funct,
    DOUBLE        **deriv,
    DOUBLE        **derxy,
    DOUBLE        **vderxy,
    DOUBLE        **xjm,
    DOUBLE        **wa1,
    DOUBLE        **xyze,
    DOUBLE        **sigmaint
    )
{
  INT     i,j,node,lr,ls,lt;   /* some counters */
  INT     iel,nir,nis,nit;     /* number of nodes/integr. points */
  INT     intc,icode;          /* flags */
  INT     actmat;              /* actual material number */
  INT     ntyp;                /* flag for element type */
  INT     iv;                  /* counter for GAUSS points */
  DOUBLE  preint,det,val;      /* element values */
  DOUBLE  e1,e2,e3,r,s,t;         
  DOUBLE  xgr[MAXGAUSS];       /* local r coords of gauss points */
  DOUBLE  xgs[MAXGAUSS];       /* local s coords of gauss points */
  DOUBLE  xgt[MAXGAUSS];       /* local t coords of gauss points */
  DIS_TYP typ;	               /* element displacement type  */
  DOUBLE  dens,visc,twovisc;   /* material parameters */
  DOUBLE  fpar[MAXGAUSS];      /* working array */
  NODE   *actnode;             /* actual node */


#ifdef DEBUG 
  dstrc_enter("f3_calelestress");
#endif


  /* initialisation */
  iel=ele->numnp;
  actmat=ele->mat-1;
  dens = mat[actmat].m.fluid->density;
  visc = mat[actmat].m.fluid->viscosity*dens; /* here we need dynamic viscosity! */
  ntyp = ele->e.f3->ntyp; 
  typ  = ele->distyp;


  switch(viscstr)
  {
    case 0: /* only real pressure */
      for (i=0;i<iel;i++)
      {
        actnode=ele->node[i];
        ele->e.f3->stress_ND.a.da[i][0]=-actnode->sol_increment.a.da[3][3]*dens;
        ele->e.f3->stress_ND.a.da[i][1]=-actnode->sol_increment.a.da[3][3]*dens;
        ele->e.f3->stress_ND.a.da[i][2]=-actnode->sol_increment.a.da[3][3]*dens;
        ele->e.f3->stress_ND.a.da[i][3]= ZERO;
        ele->e.f3->stress_ND.a.da[i][4]= ZERO;
        ele->e.f3->stress_ND.a.da[i][5]= ZERO;
      }
      break;


    case 1: /* real pressure + viscose stresses */
      /* sigma = -p_real*I + 2*nue * eps(u) */ 
      /* calculate deformation velocity tensor */
      /*------------------------------------------------ get integraton data  */
      switch (ntyp)
      {
        case 1:  /* hex - element */
          icode = 2;
          nir = ele->e.f3->nGP[0];
          nis = ele->e.f3->nGP[1];
          nit = ele->e.f3->nGP[2];
          break;

        case 2: /* tet - element */  
          if (iel>4)
            icode   = 2;
          /* initialise integration */
          nir  = ele->e.f3->nGP[0];
          nis  = 1;
          nit  = 1; 
          intc = ele->e.f3->nGP[1];  
          break;
        default:
          dserror("ntyp unknown!");
      } /* end switch(ntyp) */


      /* set element velocities, real pressure and coordinates */
      for (j=0;j<iel;j++)
      {
        evel[0][j] = ele->node[j]->sol_increment.a.da[3][0];
        evel[1][j] = ele->node[j]->sol_increment.a.da[3][1];
        evel[2][j] = ele->node[j]->sol_increment.a.da[3][2];
        epre[j]    = ele->node[j]->sol_increment.a.da[3][3]*dens;
        xyze[0][j] = ele->node[j]->x[0];
        xyze[1][j] = ele->node[j]->x[1];
        xyze[2][j] = ele->node[j]->x[2];
      }/*end for (j=0;j<iel;j++) */


      /* loop over integration points */
      iv=0;
      twovisc=TWO*visc;
      for (lr=0;lr<nir;lr++)
      {    
        for (ls=0;ls<nis;ls++)
        {
          for (lt=0;lt<nit;lt++)
          {

            /* get values of  shape functions and their derivatives */
            switch(ntyp)  
            {
              case 1:   /* --> hex - element */
                e1   = data->qxg[lr][nir-1];
                e2   = data->qxg[ls][nis-1];
                e3   = data->qxg[lt][nit-1];
                f3_hex(funct,deriv,NULL,e1,e2,e3,typ,icode);
                xgr[lr] = e1;
                xgs[ls] = e2;
                xgt[lt] = e3;
                break;
              case 2:   /* --> tet - element */		  	
                e1   = data->txgr[lr][intc];
                e2   = data->txgs[lr][intc];
                e3   = data->txgt[lr][intc]; 
                f3_tet(funct,deriv,NULL,e1,e2,e3,typ,icode); 
                break;
              default:
                dserror("ntyp unknown!");
            } /* end switch (ntyp) */


            /* compute Jacobian matrix */
            f3_jaco(funct,deriv,xjm,&det,ele,iel);

            /* compute global derivates */
            f3_gder(derxy,deriv,xjm,wa1,det,iel);

            /* get velocity  derivatives at integration point */
            f3_vder(vderxy,derxy,evel,iel);

            /* get pressure at integration point */         
            f3_prei(&preint,funct,epre,iel);


            /*
                        | Ux,x    Ux,y    Ux,z |
                        |                      |
               vderxy = | Uy,x    Uy,y    Uy,z |
                        |                      |
                        | Uz,x    Uz,y    Uz,z |


                            | Ux,x+Ux,x   Ux,y+Uy,x   Ux,z+Uz,x |
                            |                                   |
               eps(u) = 1/2 |             Uy,y+Uy,y   Uy,z+Uz,y |
                            |                                   |
                            |   symm.                 Uz,z+Uz,z |


               SIGMA = -p_real*I + 2*nue * eps(u)
            */


            sigmaint[iv][0] = -preint + twovisc * vderxy[0][0];
            sigmaint[iv][1] = -preint + twovisc * vderxy[1][1];
            sigmaint[iv][2] = -preint + twovisc * vderxy[2][2];
            sigmaint[iv][3] = (vderxy[0][1] + vderxy[1][0])*visc;
            sigmaint[iv][4] = (vderxy[1][2] + vderxy[2][1])*visc;
            sigmaint[iv][5] = (vderxy[0][2] + vderxy[2][0])*visc;

            iv++;			     

          } /* end loop over nit */
        } /* end loop over nis */
      } /* end loop over nir */


      /* extrapolate stresses to the nodes */
      f3_sext(
          ele->e.f3->stress_ND.a.da,
          sigmaint,
          xgr,
          xgs,
          xgt,
          nir,
          nis,
          nit,
          iel);
      break;


    default:
      dserror("Wrong flag for viscstr !!");
      break;

  } /* end of switch(viscstr) */

#ifdef DEBUG 
  dstrc_exit();
#endif

  return; 
} /* end of f3_calelestress */



/*!----------------------------------------------------------------------
  \brief returns local coordinates of element nodes in rst-coordinates

  <pre>                                                          mn 03/04
  This routine returns local coordinates of element nodes in rst-coordinates 
  for a 3D fluid element.

  </pre>
  \param   node     INT  (i)   element node
  \param    irs     INT  (i)   flag for r,s or t

  \warning There is nothing special to this routine

  \return DOUBLE local coordinates of element node

  \sa calling: ---; called by: ---

 *----------------------------------------------------------------------*/
DOUBLE f3_rsn (
    INT   node,
    INT   irs
    )
{

  static DOUBLE  xh8[24] = { 
     1., 1.,-1.,-1., 1., 1.,-1.,-1.,
    -1., 1., 1.,-1.,-1., 1., 1.,-1.,
    -1.,-1.,-1.,-1., 1., 1., 1., 1.
  };

  static DOUBLE xh20[36] = { 
     1.,  0., -1.,  0.,  1.,  0., -1.,  0.,  1.,  1., -1., -1., 
     0.,  1.,  0., -1.,  0.,  1.,  0., -1., -1.,  1.,  1., -1., 
    -1., -1., -1., -1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.,  0. 
  };

  INT inode;
  DOUBLE ret_val;


#ifdef DEBUG 
  dstrc_enter("f3_rsn");
#endif


  if(node<=8)
  {
    inode = node-1 + (irs-1)*8;
    ret_val=xh8[inode];
  }

  else if(node<=20)
  {
    inode=((node-8)-1)+ (irs-1)*12;
    ret_val=xh20[inode];
  }

  else
  {
    dserror("unknown number of nodes in hex element");
  }
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return(ret_val);
} /* end of f3_rsn */




/*!----------------------------------------------------------------------
  \brief calculate nth order legendre polynomial of degree 'n' at 'z'

  <pre>                                                          mn 03/04
  This routine calculates nth order legendre polynomial of degree 'n' at 'z'
  for a 3D fluid element.

  </pre>
  \param       i     INT  (i)
  \param       n     INT  (i) order legendre polynomial
  \param      zr     INT* (i)
  \param       z  DOUBLE  (i) z-coordinate
  \param   value  DOUBLE* (o) value at z

  \warning There is nothing special to this routine

  \return void                                               

  \sa calling: ---; called by: --- 

 *----------------------------------------------------------------------*/
void f3_lgpl (
    INT         i,
    INT         n,
    DOUBLE    *zr,  
    DOUBLE      z,
    DOUBLE *value
    )
{
  /*----------------------------------------------------------------------*/
  INT j;
  DOUBLE zi, zj;
#ifdef DEBUG 
  dstrc_enter("f3_lgpl");
#endif
  /*----------------------------------------------------------------------*/
  zi = zr[i];
  *value = 1.;
  for (j=0; j<n; j++)
  {
    zj = zr[j];
    if(j==i) continue;
    *value *= (z-zj)/(zi-zj);
  }
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of f3_lgpl */





/*!----------------------------------------------------------------------
  \brief subroutine f3_hxsm

  <pre>                                                          mn 03/04
  subroutine f3_hxsm

  </pre>

  \param   nir,nis,nit     INT        (i) num GP in r/s/t direction
  \param      rk,sk,tk     DOUBLE     (i) r,s,t -coordinates
  \param             f     DOUBLE[][] (i) original values on g.p.
  \param            fp     DOUBLE*    (o) extrapolated values
  \param   xgr,xgs,xgt     DOUBLE*    (i) coordinates of g.p.

  \warning There is nothing special to this routine

  \return void                                               

  \sa calling: ---; called by: ---

 *----------------------------------------------------------------------*/
void f3_hxsm (
    INT nir,
    INT nis,
    INT nit,  
    DOUBLE rk,
    DOUBLE sk,
    DOUBLE tk,
    DOUBLE f[6][27],
    DOUBLE *fp,
    DOUBLE *xgr,
    DOUBLE *xgs,
    DOUBLE *xgt
    )
{

  INT i, j, k, ns, kkk, ngp;
  DOUBLE xlr, xls, xlt;

#ifdef DEBUG 
  dstrc_enter("f3_hxsm");
#endif

  kkk=6;
  /*----------------------------------------------------------------------*/
  for (i=0; i<6; i++) fp[i] = 0.0;
  /*----------------------------------------------------------------------*/
  ngp=0;
  /*----------------------------------------------------------------------*/
  for (i=0; i<nir; i++) {
    f3_lgpl (i,nir,xgr,rk,&xlr);
    for (j=0; j<nis; j++) {
      f3_lgpl (j,nis,xgs,sk,&xls);
      for (k=0; k<nit; k++) {
        f3_lgpl (k,nit,xgt,tk,&xlt);
        for (ns=0; ns<kkk; ns++)  fp[ns] += xlr*xls*xlt*f[ns][ngp];
        ngp = ngp + 1;
      }
    }
  }

#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of f3_hxsm */




/*!----------------------------------------------------------------------
  \brief extrapolation of stress from gauss points to nodal points

  <pre>                                                          mn 03/04
  This routine extrapolates stresses from gauss points to nodal points 
  for a 3D fluid element.


  </pre>
  \param     nostrs  DOUBLE**  (o) element stresses extrapolated to the nodes
  \param   gpstress  DOUBLE**  (i) element stresses of integration points
  \param        xgr  DOUBLE*   (i) local rst-coordinates of integration points
  \param        xgs  DOUBLE*   (i) ..
  \param        xgt  DOUBLE*   (i) ..
  \param        nir     INT    (i) number of integration points in rst-direction
  \param        nis     INT    (i) ..
  \param        nit     INT    (i) ..
  \param        iel     INT    (i) number of element nodes

  \warning There is nothing special to this routine
  \return void                                               
  \sa calling: ---; called by: ---

 *----------------------------------------------------------------------*/
void f3_sext(
    DOUBLE  **nostrs,
    DOUBLE  **gpstress,
    DOUBLE   *xgr,
    DOUBLE   *xgs,
    DOUBLE   *xgt,
    INT       nir,
    INT       nis, 
    INT       nit,
    INT       iel
    )
{

  INT        nn,i,j,ngp;
  DOUBLE     cnp1, cnp2, cnp3, det;
  /* gausspoint stresses */
  DOUBLE     fgp[6][27];


#ifdef DEBUG 
  dstrc_enter("f3_sext");
#endif


  ngp = nir*nis*nit;

  for (i=0; i<6; i++)
  {
    for (j=0; j<27; j++)
    {
      fgp[i][j] = 0.0;
    }
  }

  for (i=0; i<ngp; i++)
  {
    for (j=0; j<6; j++)
    {
      fgp[j][i] = gpstress[i][j];
    }
  }


  /*----------------------------------------------------------------------*/
  for (nn=1; nn<=iel; nn++)
  {
    cnp1 = f3_rsn (nn,1);
    cnp2 = f3_rsn (nn,2);
    cnp3 = f3_rsn (nn,3);

    f3_hxsm (nir,nis,nit,cnp1,cnp2,cnp3,fgp,nostrs[nn-1],&(xgr[0]),&(xgs[0]),&(xgt[0]));
  }

#ifdef DEBUG 
  dstrc_exit();
#endif

  return;
} /* end of f3_sext */


#endif
/*! @} (documentation module close)*/

