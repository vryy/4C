/*!----------------------------------------------------------------------
\file
\brief AITKEN iteration

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/
/*!
\addtogroup FSI
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fsi_prototypes.h"
#ifdef D_FSI
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
extern struct _PAR   par;

static FSI_DYNAMIC    *fsidyn;

/*------------------------------------- prototypes used in here only ---*/
void fsi_omega_sg( FIELD *actfield );

/*!---------------------------------------------------------------------*
\brief compute relaxation parameter via steepest descent method

<pre>                                                            cf 08/03

line numbers / *=== 5.o ===* / following MOK's Algorithm 6.4

see Dissertation of D.P.MOK, p. 125 ff

</pre>

\param *structfield   FIELD	     (i)   structural field
\param *fsidyn 	      FSI_DYNAMIC    (i)
\param  itnum         INT            (i)   actual iteration step
\return void

------------------------------------------------------------------------*/
void fsi_gradient(
                  FIELD          *alefield,
                  FIELD          *structfield,
                  FIELD          *fluidfield,
		  INT             numfa,
		  INT             numff,
		  INT             numfs
	         )
{
INT       i,j;           /* counters                                    */
INT       numveldof;     /* number of velocity dofs                     */
INT       numnp_fluid;   /* number of fluid nodes                       */

DOUBLE    actr;          /* (artificial) mesh displacement              */
DOUBLE    initr;         /* initial mesh position                       */
DOUBLE    dt;            /* time step size                              */

NODE     *actfnode;      /* actual fluid node                           */
GNODE    *actfgnode;     /* actual fluid gnode                          */
NODE     *actanode;      /* actual ale node                             */
ALE_DYNAMIC    *adyn;
FLUID_DYNAMIC  *fdyn;
STRUCT_DYNAMIC *sdyn;

#ifdef DEBUG
dstrc_enter("fsi_gradient");
#endif

/*--------------------------------------------- set dynamics pointer ---*/
adyn   = alldyn[genprob.numaf].adyn;
fdyn   = alldyn[genprob.numff].fdyn;
sdyn   = alldyn[genprob.numsf].sdyn;
fsidyn = alldyn[3].fsidyn;

numnp_fluid  = fluidfield->dis[0].numnp;
numveldof    = fdyn->numdf - 1;
dt           = fsidyn->dt;

/*---------------------------------------------------- screen report ---*/
if (par.myrank==0)
   printf("Performing Steepest Descent Method for relaxation parameter\n");

/*=== 5.a ===*/
/*-- calculate residuum g_i and write it to sol_mf[6][i] of structfield */
solserv_sol_zero(structfield,0,node_array_sol_mf,6);
solserv_sol_add(structfield, 0, node_array_sol_mf, node_array_sol_mf, 0, 6, 1.0);
solserv_sol_add(structfield, 0, node_array_sol_mf, node_array_sol_mf, 1, 6,-1.0);

/*=== 5.b ===*/
/*--------------------- solve Ale mesh with sol_mf[6][i] used as DBC ---*/
fsi_ale(alefield, 6);
/* note: The ale solution of the auxiliary problem is always performed
         linear. */

/*=== 5.c,d ===*/
/*------------------------------------- create some additional space ---*/
solserv_sol_zero(fluidfield,0,node_array_sol_increment,7);
solserv_sol_zero(alefield,0,node_array_sol_mf,2);
/*------------------------------------- mesh velocity to fluid field ---*/
for (i=0;i<numnp_fluid;i++) /*----------------------- loop all nodes ---*/
{
   actfnode  = &(fluidfield->dis[0].node[i]);
   actfgnode = actfnode->gnode;
   actanode  = actfgnode->mfcpnode[numfa];

   if (actanode==NULL) continue;
   for (j=0;j<numveldof;j++)
   {
       actr   = actanode->sol_increment.a.da[0][j];
       initr  = actanode->x[j];
       /* grid velocity */
       actfnode->sol_increment.a.da[4][j] = actr/dt;
       /* grid position */
       actanode->sol_mf.a.da[2][j] = actr + initr;
   } /* end of loop over vel dofs */
} /* end of loop over all nodes */

/*=== 5.e,f,g ===*/
/*------------------------------------------------------ solve fluid ---*/
fsi_fluid(fluidfield,6);

/*=== 5.h ===*/
/*-------------------------------------------------- solve structure ---*/
fsi_struct(structfield,6,0);

/*=== 5.i ===*/
/*----------------------------------- determine relaxation parameter ---*/
fsi_omega_sg(structfield);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fsi_gradient */


/*!---------------------------------------------------------------------*
\brief evaluate relaxation parameter in steepest descent method

<pre>                                                            cf 08/03

this routine performes the final evaluation of the relaxation parameter
using the steepest descent method.

omega = (g_i^T * g_i) / (g_i^T * (-d_Gamma + g_i) )

g_i 		stored on node->sol_mf.a.da[6][i]
d_Gamma		stored on node->sol.a.da[8][i]

NOTE: This routine is called by 'fsi_gradient()' above only, where the
      static pointer fsidyn is set appropriately.

</pre>

\param *actfield   FIELD	     (i)   structural field
\return void

------------------------------------------------------------------------*/
void fsi_omega_sg( FIELD *actfield )
{
INT       i,j;           /* counters                                    */
INT       numdf;         /* degrees of freedom per node                 */
INT       numnp_total;   /* number of nodes in actfield                 */

DOUBLE    prod1;         /* auxiliary value, scalar product             */
DOUBLE    prod2;         /* auxiliary value, scalar product             */

NODE     *actnode;       /* actual structure node                       */
NODE     *actfnode;      /* actual fluid node                           */
GNODE    *actgnode;      /* actual structure geometry node              */
GNODE    *actfgnode;     /* actual fluid gnode                          */

#ifdef DEBUG
dstrc_enter("fsi_gradient");
#endif

/*----------------------------------------------------------------------*/
numnp_total  = actfield->dis[0].numnp;

prod1 = 0.0;
prod2 = 0.0;
/*------------------------- loop all nodes and look for coupling dbc ---*/
for (i=0;i<numnp_total;i++)
{
   actnode  = &(actfield->dis[0].node[i]);
   actgnode = actnode->gnode;
   actfnode = actgnode->mfcpnode[1];
   actfgnode = actfnode->gnode;
   if (actfgnode->dirich==NULL)  /* no dirichlet boundary condition */
         continue;
   switch(actfgnode->dirich->dirich_type)
   {
   case dirich_none:            /* 'ordinary' dbc */
         continue;
   case dirich_FSI:             /* coupling degree of freedom */
      numdf = actnode->numdf;
      for (j=0;j<numdf;j++)     /* loop dofs */
      {
         prod1 += actnode->sol_mf.a.da[6][j] * actnode->sol_mf.a.da[6][j];
	 prod2 += actnode->sol_mf.a.da[6][j] * ( actnode->sol_mf.a.da[6][j]
	                                       - actnode->sol.a.da[8][j] );
      }
   break;
   default:
      dserror("dirch_type unknown!\n");
   } /* end switch */
} /*end loop over nodes */

if( prod2 != 0.0 )
   fsidyn->relax = prod1/prod2;
else
   fsidyn->relax = 100.0;

/*------------------------------------------------- output to the screen */
if (par.myrank==0)
{
printf("          RELAX = %.5lf\n",fsidyn->relax);
printf("\n");
}

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fsi_omega_sg */


#endif

/*! @} (documentation module close)*/

