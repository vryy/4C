/*!----------------------------------------------------------------------
\file
\brief chimera_boundary_update.c

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_CHIMERA

/*!
\addtogroup CHIMERA
*//*! @{ (documentation module open)*/



#include "../headers/standardtypes.h"
#include "chimera.h"






/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
  *----------------------------------------------------------------------*/
extern struct _GENPROB               genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD                *field;
extern  struct _CHIMERA_DATA        *chm_data;





/*!----------------------------------------------------------------------
\brief perform boundary data exchange between discretizations

<pre>                                                            irhan 09/04

Zu gegebenem Chimera-Dirichletknoten wird der
Wert aus dem Hintergrundgitter interpoliert und zugewiesen.

Es k"onnen Dreicksgitter und Rechtecksgitter als Hintergrund
verwendet werden.

Schritte:
1) Brute-force-Suche nach parent-Element
2) Interpolation der aktuellen Dirichletrandbedingung
*/
/*
Peter Gamnitzer
last change: 2004-02-10

The value at a given Chimera-Dirichlet-node is interpolated
from the background subdomain,

quad4- and tri3-elements are possible.

steps:
1) Brute-force-search for parent-element
2) Interpolation of the Dirichlet-values

</pre>

*----------------------------------------------------------------------*/
DOUBLE chimera_boundary_update(
  NODE      *actnode,
  GNODE     *actgnode,
  INT        pos,
  INT        n_backgrounddis
  )
{
  INT          j,rr;
  DISCRET     *backgrddis;
  DOUBLE       x,y,xi,zeta,temp,tempx,tempy,Nenner;
  DOUBLE       N0[2],N1[2],N2[2],N3[2];
  DOUBLE       u0[3],u1[3],u2[3],u3[3];
  DOUBLE       steady;
  ELEMENT     *parent;
  
#ifdef DEBUG
  dstrc_enter("chimera_boundary_update");
#endif
/*----------------------------------------------------------------------*/

    backgrddis = &(field[genprob.numff].dis[n_backgrounddis]);
    x = actnode->x[0];
    y = actnode->x[1];

    /* Rechteckselemente */
    if (field[genprob.numff].dis[n_backgrounddis].element[0].distyp == quad4)
    {
      /* Brute-force Suche nach parent */
      parent = chimera_search(x,y,backgrddis,n_backgrounddis,chm_data[0].search_TOL,chm_data[0].search_action);

      /* Funktionswerte an den Ecken */
      for (rr=0; rr<3; rr++)
      {
        u0[rr] = parent->node[0]->sol_increment.a.da[pos][rr];
        u1[rr] = parent->node[1]->sol_increment.a.da[pos][rr];
        u2[rr] = parent->node[2]->sol_increment.a.da[pos][rr];
        u3[rr] = parent->node[3]->sol_increment.a.da[pos][rr];
      }

      /* Interpolation der aktuellen Randwerte */
      for (rr=0; rr<2; rr++)
      {
        N0[rr] = (parent->node[0])->x[rr];
        N1[rr] = (parent->node[1])->x[rr];
        N3[rr] = (parent->node[3])->x[rr];
      }
      tempx = N1[0]*N1[0]+N0[0]*N0[0]+N1[1]*N1[1]+N0[1]*N0[1]-2*(N0[1]*N1[1]+N0[0]*N1[0]);
      xi = ((x-N0[0])*(N1[0]-N0[0])+(y-N0[1])*(N1[1]-N0[1]))/tempx;
      tempy = N3[0]*N3[0]+N0[0]*N0[0]+N3[1]*N3[1]+N0[1]*N0[1]-2*(N0[1]*N3[1]+N0[0]*N3[0]);
      zeta = ((x-N0[0])*(N3[0]-N0[0])+(y-N0[1])*(N3[1]-N0[1]))/tempy;
      /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
      for (j=0; j<actnode->numdf-1; j++) /* loop dofs */
      {
        actnode->sol_increment.a.da[pos][j] = (u0[j]+(u1[j]-u0[j])*xi+(u3[j]-u0[j])*zeta +
                                               (u2[j]+u0[j]-u1[j]-u3[j])*xi*zeta);
/*	    printf("1-xi-zeta*(1-xi) %e  zeta*1-xi*zeta %e  xi*(1-zeta) %e xi %e  zeta %e\n",1-xi-zeta*(1-xi), zeta-xi*zeta,xi*(1-zeta),xi,zeta);*/
/*	    actnode->sol_increment.a.da[pos+2][j] =u0[j]*(tempx*tempy-xi*tempy-zeta*(tempx-xi))+u1[j]*xi*(tempy-zeta)+u2[j]*xi*zeta+u3[j]*zeta*(tempx-xi);
	    actnode->sol_increment.a.da[pos+2][j] =actnode->sol_increment.a.da[pos+2][j]/(tempx*tempy);*/
/*	    actnode->sol_increment.a.da[pos+2][j] =u0[j]*(1-xi-zeta*(1-xi))+u1[j]*xi*(1-zeta)+u2[j]*xi*zeta+u3[j]*zeta*(1-xi);*/

/*	    printf("xi  %e   zeta %e    actnode->sol_increment.a.da[pos+2][j]  %e\n",xi,zeta, actnode->sol_increment.a.da[pos+2][j]);*/
	}
      if (actnode[0].gnode[0].chi_pres_coupling_point == chimera_yes)
      {
        j = 2;
        actnode->sol_increment.a.da[pos][j] = 1*(u0[j]+(u1[j]-u0[j])*xi+(u3[j]-u0[j])*zeta +
                                                 (u2[j]+u0[j]-u1[j]-u3[j])*xi*zeta);
      }
      /* if (n_backgrounddis==0) {
	    printf("Identifikation dis %d\n",backgrddis[0].numele);
	    printf("Hintergrundupdate\n");
	    printf("%5f   %5f   %5f\n",actnode->sol_increment.a.da[pos+2][0],actnode->sol_increment.a.da[pos+2][1],actnode->sol_increment.a.da[pos+2][2]);
	    }*/
    }

    /* Dreieckselemente */
    if (field[genprob.numff].dis[n_backgrounddis].element[0].distyp == tri3)
    {
      /* brute-force Suche nach parent */
      parent = chimera_search(x,y,backgrddis,n_backgrounddis,chm_data[0].search_TOL,chm_data[0].search_action);
      /* Funktionswerte an den Ecken */
      for (rr=0; rr<3; rr++)
      {
        u0[rr] = parent->node[0]->sol_increment.a.da[pos][rr];
        u1[rr] = parent->node[1]->sol_increment.a.da[pos][rr];
        u2[rr] = parent->node[2]->sol_increment.a.da[pos][rr];
      }

      /* Interpolation der aktuellen Randwerte */
      for (rr=0; rr<2; rr++)
      {
        N0[rr] = (parent->node[0])->x[rr];
        N1[rr] = (parent->node[1])->x[rr];
        N2[rr] = (parent->node[2])->x[rr];
      }
      temp = 1/((N1[0]-N0[0])*(N1[1]-N2[1])-(N1[0]-N2[0])*(N1[1]-N0[1]));
      xi = (x*(N0[1]-N2[1])+y*(N2[0]-N0[0])+N0[0]*N2[1]-N2[0]*N0[1])*temp;
      temp = 1/((N2[0]-N0[0])*(N2[1]-N1[1])-(N2[0]-N1[0])*(N2[1]-N0[1]));
      zeta = (x*(N0[1]-N1[1])+y*(N1[0]-N0[0])+N0[0]*N1[1]-N1[0]*N0[1])*temp;
      for (j=0; j<actnode->numdf-1; j++) /* loop dofs*/
      {
        actnode->sol_increment.a.da[pos][j] = u1[j]*xi+u2[j]*zeta+u0[j]*(1-xi-zeta);
      }
      if (actnode[0].gnode[0].chi_pres_coupling_point == chimera_yes)
      {
        j = 2;
        actnode->sol_increment.a.da[pos][j] = u1[j]*xi+u2[j]*zeta+u0[j]*(1-xi-zeta);
      }
/*	printf("%5f   %5f   %5f\n",actnode->sol_increment.a.da[pos+2][0],actnode->sol_increment.a.da[pos+2][1],actnode->sol_increment.a.da[pos+2][2]);*/
    }

    /* Konvergenz"uberpr"ufung */
    /* convergence-check (used as a stopping criteria for the subdomain iteration) */
    /* "relativer" Zuwachs */
    steady = 0;
    temp = 0;
    Nenner = 0;
    for (j=0; j<actnode->numdf-1; j++)
    {
      if(actnode->sol_increment.a.da[pos][j] != 0)
      {
        temp = (actnode->sol_increment.a.da[pos][j]-actnode->sol_increment.a.da[pos+1][j]);
        temp = temp*temp;
        Nenner = Nenner + (actnode->sol_increment.a.da[pos+1][j]*actnode->sol_increment.a.da[pos+1][j]);
        steady = steady + temp;
      }
    }
    if(Nenner != 0)
    {
      steady = steady/Nenner;
    }
    steady = sqrt(steady);
    /*    printf("chm_data[0].rel_err[%d]  %5f       steady %5f\n",1-n_backgrounddis,chm_data[0].rel_err[1-n_backgrounddis],steady);*/
    if (steady > chm_data[0].rel_err[1-n_backgrounddis])
    {
      chm_data[0].rel_err[1-n_backgrounddis] = steady;
    }
    /* absoluter Zuwachs*/
    steady = 0;
    temp = 0;
    for (j=0; j<actnode->numdf-1; j++)
    {
      temp = (actnode->sol_increment.a.da[pos][j]-actnode->sol_increment.a.da[pos+1][j]);
      temp = temp*temp;
      steady = steady + temp;
    }
    steady = sqrt(steady);
    /*printf("chm_data[0].abs_err[%d]  %5f       steady %5f\n",1-n_backgrounddis,chm_data[0].abs_err[1-n_backgrounddis],steady);*/
    if (steady>chm_data[0].abs_err[1-n_backgrounddis])
    {
      chm_data[0].abs_err[1-n_backgrounddis] = steady;
    }

/*----------------------------------------------------------------------*/    
#ifdef DEBUG
    dstrc_exit();
#endif
    
    return(steady);
} /* end of chimera_boundary_update */
/*! @} (documentation module close)*/
#endif
