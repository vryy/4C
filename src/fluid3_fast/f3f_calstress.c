/*-----------------------------------------------------------------------*/
/*!
\file
\brief calculation of fluid3_fast stresses

  Very detailed description.

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/

/*!
\addtogroup Fluid3_fast
*//*! @{ (documentation module open)*/


#ifdef D_FLUID3_F


#include "../headers/standardtypes.h"
#include "f3f_prototypes.h"
#include "../fluid3/fluid3.h"
/*----------------------------------------------------------------------*
  |                                                     m.gee 06/01    |
  | vector of material laws                                            |
  | defined in global_control.c                                        |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | pointer to allocate dynamic variables if needed                      |
  | dedfined in global_control.c                                         |
  | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;

/*!----------------------------------------------------------------------
\brief positions of physical values in node arrays

<pre>                                                        chfoe 11/04

This structure contains the positions of the various fluid solutions 
within the nodal array of sol_increment.a.da[ipos][dim].

extern variable defined in fluid_service.c
</pre>

------------------------------------------------------------------------*/
extern struct _FLUID_POSITION ipos;



/*-----------------------------------------------------------------------*/
/*!
  \brief calculation of fluid3_fast stresses

  calculation of the stress tensor sigma, which es stored at the nodes

  sigma = -p_real*I + 2*nue * eps(u)

  \param viscstr         INT      (i) include viscose stresses yes/no
  \param ele[LOOPL]      ELEMENT  (i) the set of elements
  \param evel           *DOUBLE   (o) vels
  \param epre           *DOUBLE   (o) pres
  \param funct          *DOUBLE   (o) shape functions
  \param deriv          *DOUBLE   (o) deriv. of shape funcs
  \param derxy          *DOUBLE   (o) global derivatives
  \param xjm            *DOUBLE   (o) jacobian matrix
  \param wa1            *DOUBLE   (o) working array
  \param elecord        *DOUBLE   (o) vector containing nodal coordinates
  \param sigint         *DOUBLE   (o) stresses at GAUSS point
  \param nostr          *DOUBLE   (o) stresses at nodes
  \param sizevec[6]      INT      (i) some sizes

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3fcalelestress(
    INT             viscstr,
    ELEMENT        *ele[LOOPL],
    DOUBLE         *evel,
    DOUBLE         *epre,
    DOUBLE         *funct,
    DOUBLE         *deriv,
    DOUBLE         *derxy,
    DOUBLE         *vderxy,
    DOUBLE         *xjm,
    DOUBLE         *wa1,
    DOUBLE         *elecord,
    DOUBLE         *sigint,
    DOUBLE         *nostr,
    INT             sizevec[6]
    )
{

  INT     i,j,l,lr,ls,lt,nn;   /* some counters */
  INT     iel,nir=0,nis=0,nit=0;     /* number of nodes/integr. points */
  INT     intc=0,inttyp;                /* flags */
  INT     icode = 0;           /* flags */
  INT     actmat;              /* actual material number */
  INT     iv;                  /* counter for GAUSS points */
  INT     velnp;               /* position flag for sol_increment field */
  DOUBLE  preint[LOOPL],det[LOOPL];          /* element values */
  DOUBLE  e1,e2,e3;
  DOUBLE  xgr[MAXGAUSS];       /* local r coords of gauss points */
  DOUBLE  xgs[MAXGAUSS];       /* local s coords of gauss points */
  DOUBLE  xgt[MAXGAUSS];       /* local t coords of gauss points */
  DIS_TYP typ;                 /* element displacement type  */
  DOUBLE  dens,visc,twovisc;   /* material parameters */
  NODE   *actnode;             /* actual node */

  DOUBLE  **stressND;

#ifdef D_FSI
  NODE   *actfnode;             /* actual node */
  NODE   *actanode;             /* actual node */
  GNODE  *actfgnode;             /* actual node */
#endif

  FLUID_DATA     *data;
  FLUID_DYNAMIC  *fdyn;

#ifdef DEBUG
  dstrc_enter("f3fcalelestress");
#endif


  /* initialisation */
  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;
  iel    = ele[0]->numnp;
  actmat = ele[0]->mat-1;
  dens   = mat[actmat].m.fluid->density;
  visc   = mat[actmat].m.fluid->viscosity*dens; /* here we need dynamic viscosity! */
  typ    = ele[0]->distyp;

  velnp = ipos.velnp;

  switch(viscstr)
  {
    case 0:
      /* only real pressure */
      for (l=0;l<sizevec[4];l++)
      {
        for (i=0;i<iel;i++)
        {
          actnode=ele[l]->node[i];
          ele[l]->e.f3->stress_ND.a.da[i][0]=-actnode->sol_increment.a.da[velnp][3]*dens;
          ele[l]->e.f3->stress_ND.a.da[i][1]=-actnode->sol_increment.a.da[velnp][3]*dens;
          ele[l]->e.f3->stress_ND.a.da[i][2]=-actnode->sol_increment.a.da[velnp][3]*dens;
          ele[l]->e.f3->stress_ND.a.da[i][3]= ZERO;
          ele[l]->e.f3->stress_ND.a.da[i][4]= ZERO;
          ele[l]->e.f3->stress_ND.a.da[i][5]= ZERO;
        }
      }
      break;


    case 1: /* real pressure + viscose stresses */

      /* sigma = -p_real*I + 2*nue * eps(u) */

      /* get integraton data  */
      switch (typ)
      {
        case hex8: case hex20: case hex27:   /* --> hex - element */
          icode = 2;
          nir = ele[0]->e.f3->nGP[0];
          nis = ele[0]->e.f3->nGP[1];
          nit = ele[0]->e.f3->nGP[2];
          intc = 0;
          break;

        case tet4: case tet10:   /* --> tet - element */
          if (iel>4)
            icode   = 2;
          /* initialise integration */
          nir  = ele[0]->e.f3->nGP[0];
          nis  = 1;
          nit  = 1;
          intc = ele[0]->e.f3->nGP[1];
          break;

        default:
          dserror("typ unknown!");
      } /* end switch(typ) */


      /* set element velocities, real pressure and coordinates */
      for(i=0;i<sizevec[1];i++) /* loop nodes */
      {
        for(j=0;j<sizevec[4];j++)
        {
          actnode=ele[j]->node[i];

          evel[LOOPL*i+j]                   = actnode->sol_increment.a.da[velnp][0];
          evel[sizevec[0]*LOOPL+LOOPL*i+j]  = actnode->sol_increment.a.da[velnp][1];
          evel[2*sizevec[0]*LOOPL+LOOPL*i+j]= actnode->sol_increment.a.da[velnp][2];
          epre[LOOPL*i+j]                   = actnode->sol_increment.a.da[velnp][3]*dens;
        } /*loop*/
      } /* end of loop over nodes */



      switch (ele[0]->e.f3->is_ale)
      {
        case 0:
          for (l=0; l<sizevec[1]; l++)
          {
            for(j=0;j<sizevec[4];j++)
            {
              elecord[                   LOOPL*l+j]=ele[j]->node[l]->x[0];
              elecord[  sizevec[0]*LOOPL+LOOPL*l+j]=ele[j]->node[l]->x[1];
              elecord[2*sizevec[0]*LOOPL+LOOPL*l+j]=ele[j]->node[l]->x[2];
            }
          }
          break;

#ifdef D_FSI
        case 1:
          for (l=0; l<sizevec[1]; l++)
          {
            for(j=0;j<sizevec[4];j++)
            {
              actfnode  = ele[j]->node[l];
              actfgnode = actfnode->gnode;
              actanode  = actfgnode->mfcpnode[genprob.numaf];
              elecord[                   LOOPL*l+j]
                = ele[j]->node[l]->x[0] + actanode->sol_mf.a.da[1][0];
              elecord[  sizevec[0]*LOOPL+LOOPL*l+j]
                = ele[j]->node[l]->x[1] + actanode->sol_mf.a.da[1][1];
              elecord[2*sizevec[0]*LOOPL+LOOPL*l+j]
                = ele[j]->node[l]->x[2] + actanode->sol_mf.a.da[1][2];
            }
          }
          break;
#endif

        default:
          dserror("elment flag is_ale out of range!");
      }


      /* loop over integration points */
      iv=1;
      twovisc=TWO*visc;

      for (lr=0;lr<nir;lr++)
      {
        for (ls=0;ls<nis;ls++)
        {
          for (lt=0;lt<nit;lt++)
          {

            /* get values of  shape functions and their derivatives */
            switch(typ)
            {
              case hex8: case hex20: case hex27:   /* --> hex - element */
                e1   = data->qxg[lr][nir-1];
                e2   = data->qxg[ls][nis-1];
                e3   = data->qxg[lt][nit-1];

                xgr[lr] = e1;
                xgs[ls] = e2;
                xgt[lt] = e3;

                if (typ==hex8)
                  inttyp=8;
                else if(typ==hex20)
                  inttyp=20;
                else if(typ==hex27)
                  inttyp=27;

                f3fhex(funct, deriv, NULL, &e1, &e2, &e3, &inttyp, &icode, sizevec);

                break;

              case tet4: case tet10:   /* --> tet - element */
                e1   = data->txgr[lr][intc];
                e2   = data->txgs[lr][intc];
                e3   = data->txgt[lr][intc];

                if (typ==tet4)
                  inttyp=4;
                else if(typ==tet10)
                  inttyp=10;

                f3ftet(funct, deriv, NULL, &e1, &e2, &e3, &inttyp, &icode, sizevec);

                break;
              default:
                dserror("typ unknown!");
            } /* end switch (typ) */


            /* compute Jacobian matrix */
            f3fjaco(funct, deriv, xjm, det, elecord, sizevec);


            /* compute global derivates */
            f3fgder(derxy, deriv, xjm, wa1, det, sizevec);

            /* get velocity (n+g,i) derivatives at integration point */
            f3fvder(vderxy, derxy, evel, sizevec);

            /* get pressure (n) at integration point */
            f3fprei(preint, funct, epre, sizevec);


            /* calculate stresses */
            f3fsigint(preint, vderxy, sigint, &visc, &iv, sizevec);

            iv++;

          } /* end loop over nit */
        } /* end loop over nis */
      } /* end loop over nir */


      /* extrapolate stresses to the nodes */
      f3fsext(
          nostr,
          sigint,
          xgr,
          xgs,
          xgt,
          nir,
          nis,
          nit,
          sizevec);

      /* copy the result back to the nodes */
      for (nn=0;nn<sizevec[1];nn++)
      {
        for (l=0;l<sizevec[4];l++)
        {
          stressND = ele[l]->e.f3->stress_ND.a.da;
          stressND[nn][0] = nostr[0*sizevec[0]*sizevec[3] + nn*sizevec[3] + l];
          stressND[nn][1] = nostr[1*sizevec[0]*sizevec[3] + nn*sizevec[3] + l];
          stressND[nn][2] = nostr[2*sizevec[0]*sizevec[3] + nn*sizevec[3] + l];
          stressND[nn][3] = nostr[3*sizevec[0]*sizevec[3] + nn*sizevec[3] + l];
          stressND[nn][4] = nostr[4*sizevec[0]*sizevec[3] + nn*sizevec[3] + l];
          stressND[nn][5] = nostr[5*sizevec[0]*sizevec[3] + nn*sizevec[3] + l];
        }
      }

      break;


    default:
      dserror("parameter viscstr out of range!\n");
  }
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3fcalelestress */




/*-----------------------------------------------------------------------*/
/*!
  \brief extrapolate stresses to the nodes

  \param nostr          *DOUBLE   (o) stresses at nodes
  \param sigint         *DOUBLE   (i) stresses at GAUSS point
  \param xgr            *DOUBLE   (i) coords of gauss points
  \param xgs            *DOUBLE   (i) coords of gauss points
  \param xgt            *DOUBLE   (i) coords of gauss points
  \param nir             INT      (i) number of gauss points
  \param nis             INT      (i) number of gauss points
  \param nit             INT      (i) number of gauss points
  \param sizevec[6]      INT      (i) some sizes

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3fsext(
    DOUBLE   *nostr,
    DOUBLE   *sigint,
    DOUBLE   *xgr,
    DOUBLE   *xgs,
    DOUBLE   *xgt,
    INT       nir,
    INT       nis,
    INT       nit,
    INT       sizevec[6]
    )

{

  INT        nn;
  DOUBLE     cnp1, cnp2, cnp3;

  INT i, j, k, l, ns, ngp;
  DOUBLE xlr, xls, xlt;

#ifdef DEBUG
  dstrc_enter("f3fsext");
#endif


  /* loop all nodes */
  for (nn=0; nn<sizevec[1]; nn++)
  {
    cnp1 = f3frsn(nn+1,1);
    cnp2 = f3frsn(nn+1,2);
    cnp3 = f3frsn(nn+1,3);


    /* initialize the nodal stresses */
    for (i=0; i<6; i++)
      for (l=0; l<sizevec[4]; l++)
        nostr[i*sizevec[0]*sizevec[3] + nn*sizevec[3] + l] = 0.0;


    ngp=0;

    /* loop all gausspoints */
    for (i=0; i<nir; i++) {
      f3flgpl(i,nir,xgr,cnp1,&xlr);

      for (j=0; j<nis; j++) {
        f3flgpl(j,nis,xgs,cnp2,&xls);

        for (k=0; k<nit; k++) {
          f3flgpl(k,nit,xgt,cnp3,&xlt);

          /* loop all stresses */
          for (ns=0; ns<6; ns++)
            for (l=0; l<sizevec[4]; l++)
              nostr[ns*sizevec[0]*sizevec[3] + nn*sizevec[3] + l]
                += xlr*xls*xlt*sigint[ns*sizevec[5]*sizevec[3] + ngp*sizevec[3] + l];

          ngp = ngp + 1;

        }
      }
    }

  } /* for (nn=0; nn<sizevec[1]; nn++) */

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3fsext */



/*-----------------------------------------------------------------------*/
/*!
  \brief returns local coordinates of element nodes in rst-coordinates

  This routine returns local coordinates of element nodes in rst-coordinates
  for a 3D fluid_fast element.

  \param node            INT      (i) element node
  \param irs             INT      (i) flag for r,s or t

  \return DOUBLE

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
DOUBLE f3frsn(
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
  dstrc_enter("f3frsn");
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
    ret_val = 0.0;
    dserror("unknown number of nodes in hex element");
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return(ret_val);
} /* end of f3frsn */




/*-----------------------------------------------------------------------*/
/*!
  \brief calculate nth order legendre polynomial of degree 'n' at 'z'

  This routine calculates nth order legendre polynomial of degree 'n' at 'z'
  for a 3D fluid element.

  \param i               INT      (i)
  \param n               INT      (i) order legendre polynomial
  \param zr              INT*     (i)
  \param z               DOUBLE   (i) z-coordinate
  \param value           DOUBLE*  (o) value at z

  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void f3flgpl(
    INT         i,
    INT         n,
    DOUBLE     *zr,
    DOUBLE      z,
    DOUBLE     *value
    )
{

  INT j;
  DOUBLE zi, zj;

#ifdef DEBUG
  dstrc_enter("f3flgpl");
#endif

  zi = zr[i];
  *value = 1.;
  for (j=0; j<n; j++)
  {
    zj = zr[j];
    if(j==i) continue;
    *value *= (z-zj)/(zi-zj);
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3flgpl */



#endif /* ifdef D_FLUID3_F */

/*! @} (documentation module close)*/


