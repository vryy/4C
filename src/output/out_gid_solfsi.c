#include "../headers/standardtypes.h"
#include "gid.h"
#include "../shell8/shell8.h"
#include "../shell9/shell9.h"
#include "../wall1/wall1.h"
#include "../brick1/brick1.h"
#include "../ale2/ale2.h"
#include "../ale3/ale3.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;   
/*----------------------------------------------------------------------*
 |                                                          mn 06/02    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
 extern struct _IO_FLAGS     ioflags;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;                      



/*----------------------------------------------------------------------*
 |  routine to write solution of a step to GID           m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_sol_fsi(FIELD *fluid, FIELD *structure)
{
INT           i,j;

static FLUID_DYNAMIC  *fdyn;
static FSI_DYNAMIC    *fsidyn;
static STRUCT_DYNAMIC *sdyn;

INTRA        *actintras;
INTRA        *actintraf;

FILE         *out = allfiles.gidres;

NODE         *actnode;
NODE         *actanode;

char         *resulttype;
char         *resultplace;
char         *rangetable;
INT           ncomponent;
char         *componentnames[18];

char          sign='"';


#ifdef DEBUG 
dstrc_enter("out_gid_sol_fsi");
#endif
/*----------------------------------------------------------------------*/
sdyn= alldyn[0].sdyn;
fdyn= alldyn[1].fdyn;
fsidyn= alldyn[3].fsidyn;
/*----------------------------------------------------------------------*/
#ifdef PARALLEL 
actintraf    = &(par.intra[genprob.numff]);
actintras    = &(par.intra[genprob.numsf]);
#else
actintraf   = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintraf) dserror("Allocation of INTRA failed");
actintraf->intra_fieldtyp = fluid;
actintraf->intra_rank     = 0;
actintraf->intra_nprocs   = 1;
actintras   = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintras) dserror("Allocation of INTRA failed");
actintras->intra_fieldtyp = structure;
actintras->intra_rank     = 0;
actintras->intra_nprocs   = 1;
#endif


if ( sdyn->step != fdyn->step)
  dserror("Something is wrong!!");

fprintf(out,"#-------------------------------------------------------------------------------\n");
fprintf(out,"# Converged Solution of timestep %d\n",fdyn->step);
fprintf(out,"#-------------------------------------------------------------------------------\n");

/* ------------- write displacements of structure field */
/*----------------------------------------------------------------------*/
if (ioflags.struct_disp_gid==1||ioflags.fluid_sol_gid==1)
{
   resulttype        = "VECTOR";
   resultplace       = "ONNODES";
   rangetable        = "standard_structure";
   ncomponent        = genprob.ndim;
   componentnames[0] = "x-displ";
   componentnames[1] = "y-displ";
   componentnames[2] = "z-displ";
   /*-------------------------------------------------------------------*/
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# RESULT DISPLACEMENTS on FIELD FSI\n");
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"RESULT %cdisplacement%c %cpcarat%c %d %s %s\n",
                                                             sign,sign,
                                                             sign,sign,
                                                             sdyn->step,
                                                             resulttype,
                                                             resultplace
                                                             );
   fprintf(out,"RESULTRANGESTABLE %c%s%c\n",
                                            sign,rangetable,sign
                                            );
   /*-------------------------------------------------------------------*/
   switch (genprob.ndim)
   {
     case 3:
       fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c,%c%s%c\n",
                                                       sign,componentnames[0],sign,
                                                       sign,componentnames[1],sign,
                                                       sign,componentnames[2],sign
                                                       );
     break;
     case 2:
       fprintf(out,"COMPONENTNAMES %c%s%c,%c%s%c\n",
                                                       sign,componentnames[0],sign,
                                                       sign,componentnames[1],sign
                                                       );
     break;
     default:
       dserror("Unknown numer of dimensions");
     break;
   }
   /*-------------------------------------------------------------------*/
   fprintf(out,"VALUES\n");
}
if (ioflags.struct_disp_gid==1)
{
   for (i=0; i<structure->dis[0].numnp; i++)
   {
      actnode = &(structure->dis[0].node[i]);
      switch (genprob.ndim)
      {
	case 3:
        fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                                                   actnode->Id+1,
                                                   actnode->sol.a.da[0][0],
                                                   actnode->sol.a.da[0][1],
                                                   actnode->sol.a.da[0][2]
                                                   );
	break;
	case 2:
        fprintf(out," %6d %18.5E %18.5E \n",
                                                   actnode->Id+1,
                                                   actnode->sol.a.da[0][0],
                                                   actnode->sol.a.da[0][1]
                                                   );
        break;
	default:
	  dserror("Unknown number of dimensions");
        break;
      }
   }


/* write displacements of fluid field */
/*-------------------------------------------------------------------*/
}
if (ioflags.fluid_sol_gid==1)
{
     for (i=0; i<fluid->dis[0].numnp; i++)
     {
       actnode = &(fluid->dis[0].node[i]);
       actanode = (fluid->dis[0].node[i].gnode->mfcpnode[genprob.numaf]);
       switch (genprob.ndim)
       {
         case 3:
           fprintf(out," %6d %18.5E %18.5E %18.5E\n",
               actnode->Id+1,
               actanode->sol.a.da[fsidyn->actpos-1][0],
               actanode->sol.a.da[fsidyn->actpos-1][1],
               actanode->sol.a.da[fsidyn->actpos-1][2]
               );
           break;
         case 2:
           fprintf(out," %6d %18.5E %18.5E \n",
               actnode->Id+1,
               actanode->sol.a.da[fsidyn->actpos-1][0],
               actanode->sol.a.da[fsidyn->actpos-1][1]
               );
           break;
         default:
           dserror("Unknown number of dimensions");
           break;
       }
     }
}
if (ioflags.struct_disp_gid==1||ioflags.fluid_sol_gid==1)
{
   fprintf(out,"END VALUES\n");
}
/*-------------------------------------------------------------------*/



/* write velocities and pressure of fluid field */
if (ioflags.fluid_sol_gid==1) 
{
  out_gid_sol("velocity",fluid,actintraf,fdyn->step,fsidyn->actpos-1);
  out_gid_sol("pressure",fluid,actintraf,fdyn->step,fsidyn->actpos-1);
}


/* write velocities, accelerations and stresses of structure field */
if (ioflags.struct_disp_gid==1)
{
  out_gid_sol("velocities",structure,actintras,sdyn->step,1);
  out_gid_sol("accelerations",structure,actintras,sdyn->step,2);
}
if (ioflags.struct_stress_gid==1)
  out_gid_sol("stress"      ,structure,actintras,sdyn->step,0);


/*----------------------------------------------------------------------*/
fflush(out);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_gid_sol_fsi */



