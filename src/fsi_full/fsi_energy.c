/*!----------------------------------------------------------------------
\file
\brief computation of fsi-interface energy

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
#ifdef D_FSI
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fsi_prototypes.h"    
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

static FSI_DYNAMIC *fsidyn;    
/*!---------------------------------------------------------------------                                         
\brief interface energy

<pre>                                                         genk 01/03

compute increment of energy transported across the fs - interface
(PIPERNO 2000)

   DESINT = [Ps(n) + Ps(n+1)]/2 * dispi 
   
   Ps(n+1) ... time discretisation of the structural interface forces
   -->  [Ps(n) + Ps(n+1)]/2 = (1-ALPHA_f)*P(n+1) + ALPHA_f*P(n)

   dispi  = sol_mf[3]    ... displacement increment (n) --> (n+1)           	      	  
   P(n+1) = sol_mf[4][j] ... coupling forces at the end of the time-step
   P(n)   = sol_mf[5][j] ... coupling forces at the beginning of the time-step
 
  DEFINT = -Pf(n+1)* [d(n+1) - d(n)]
         = -Pf(n+1)* (sol_mf[1][j] - sol_mf[2][j])
	 	 
  Pf(n+1) ... time discretisation of the fluid interface pressure
  	      Pf(n+1) = THETA*P(n+1) + (1-THETA)*P(n)
  sol_mf[1][j]  ... (relaxed) displacements of previous iteration step
  sol_mf[2][j]  ... (relaxed) displacements of previous time step
  
</pre>

\param *structfield   FIELD	     (i)   structural field
\param *fsidyn 	      FSI_DYNAMIC    (i/o)   
\param *sdyn          STRUCT_DYNAMIC (i)
\param *fdyn          FLUID_DYNAMIC  (i)
\param  init          INT            (i)  initialisation flag
\warning this function is not parallised but it's no problem!!!
\return void                                                                             

------------------------------------------------------------------------*/
void fsi_dyneint( 
                       FIELD          *structfield, 
		       INT             init
		)
{
INT             i,j;            /* some counters                        */
INT             numdf;          /* actual number of dofs                */
INT             dof;            /* actual dof                           */
static INT      numnp_total;    /* total number of structure nodes      */
static INT      numdf_total;    /* total number of structure dofs       */
static INT     *sid;            /* structural interface dofs            */
DOUBLE          desint,defint;  /* energies                             */
DOUBLE          fac;       
DOUBLE          dispi;          /* displacement increment               */
DOUBLE        **sol_mf;         /* multifield nodal solution array      */
static DOUBLE   alpha;          /* structural time integr. factor       */
static DOUBLE   oma;            /* 1-alpha                              */
static DOUBLE   theta;          /* fluid time integr. factor            */
static DOUBLE   omt;            /* 1-theta                              */
NODE           *actsnode;       /* actual structure node                */
static STRUCT_DYNAMIC *sdyn;
static FLUID_DYNAMIC  *fdyn;

#ifdef DEBUG 
dstrc_enter("fsi_dyneint");
#endif

fsidyn      = alldyn[3].fsidyn;
sdyn        = alldyn[genprob.numsf].sdyn;
fdyn        = alldyn[genprob.numff].fdyn;

/*----------------------------------------------------- set some values */
if (init==1) /* initialisation called by structure */
{
   sid         = fsidyn->sid.a.iv;
   numdf_total = fsidyn->sid.fdim;
   numnp_total = structfield->dis[0].numnp;
   alpha       = sdyn->alpha_f;
   oma         = ONE - alpha;
   desint      = ZERO;
   goto end;
}
else if (init==2) /* initialisation called by fluid */
{
   dsassert(fdyn->iop==4,
   "computation of interface energy only implemented for ONE-STEP-THETA!\n");
   theta       = fdyn->theta;
   omt         = ONE-theta;
   goto end;
}   

dsassert(fsidyn->ifsi!=3,
"computation of interface energy not implemented for fsi-algo with dt/2-shift!\n");


/*-------------------------------------------------- energy computation */
desint = ZERO;
defint = ZERO;

/*---------------------------------------------------------- loop nodes */
for (i=0;i<numnp_total;i++)
{
   actsnode  = &(structfield->dis[0].node[i]);
   numdf = actsnode->numdf; 
   sol_mf = actsnode->sol_mf.a.da;
   /*--------------------------------- loop dofs and check for coupling */
   for (j=0;j<numdf;j++)
   {
      dof = actsnode->dof[j];
      dsassert(dof<numdf_total,"dofnumber not valid!\n");
      if (sid[dof]==0) continue;
      /*----------------------------------------energy by the structure */
      dispi   = sol_mf[3][j];
      fac     = oma*sol_mf[4][j] + alpha*sol_mf[5][j];
      desint += fac*dispi;
      /*------------------------------------------- energy by the fluid */
      dispi   = sol_mf[1][j] - sol_mf[2][j];
      fac     = theta*sol_mf[4][j] + omt*sol_mf[5][j];
      defint -= fac*dispi;
   } /* end of loop over dofs */   
} /* end of loop over nodes */

/*---------------------------------- energy production at the interface */
fsidyn->deltaeint = desint + defint;

end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fsi_dyneint */

/*!---------------------------------------------------------------------                                         
\brief inteface energy check

<pre>                                                         genk 01/03

check energy production at the fs-interface
  
</pre>

\param *fsidyn 	      FSI_DYNAMIC    (i)   
\return void                                                                             

------------------------------------------------------------------------*/
void fsi_energycheck(void) 
{

#ifdef DEBUG 
dstrc_enter("fsi_energycheck");
#endif

fsidyn = alldyn[3].fsidyn;

if (fsidyn->deltaeint>fsidyn->entol)
{
   printf("\n");
   printf("##################################################\n");
   printf("SUM OF ENERGY CREATED AT FS-INTERFACE > TOLERANCE:\n");
   printf("DELTAEINT = %10.3E > ENTOL = %10.3E\n",
           fsidyn->deltaeint,fsidyn->entol);
   printf("  ==> COMPUTATION BECOMES UNSTABLE ... STOPPING \n");
   printf("##################################################\n");
   printf("\n");
   /*---------------------------------- stop calculation regularily!!! */
   fsidyn->step=fsidyn->nstep;
}

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fsi_energycheck */

#endif
/*! @} (documentation module close)*/
