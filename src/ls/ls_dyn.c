#ifdef D_LS
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "../xfem/xfem_prototypes.h"
#include "ls_prototypes.h"    



extern ALLDYNA           *alldyn;
extern struct _GENPROB    genprob;
extern struct _FIELD     *field;
extern struct _PAR        par;  
extern INT                numcurve;
extern struct _CURVE     *curve;





/************************************************************************
 ----------------------------------------- last checked by Irhan 28.04.04
 ************************************************************************/
void ls_dyn()
{
  INT                i;
  INT		     eflag;
  FIELD             *fluidfield;
  FIELD             *lsetfield;
  FLUID_DYNAMIC     *fdyn;
  LS_DYNAMIC        *lsdyn;
  INT                numfld;
  INT                numff;
  INT                numls;
  INT                actcurve;

  LS_GEN_DATA       *lsdata;
  LSFLAG             lsflag;
    
#ifdef DEBUG 
  dstrc_enter("ls_dyn");
#endif
/*----------------------------------------------------------------------*/

  /*
   * convention used at the moment =>
   * FIELD 0: levelset
   * FIELD 1: fluid
   */
  numfld=genprob.numfld;
  dsassert(numfld==2,"Two fields needed for Lset-problem!\n");
  numls       = genprob.numls;
  lsetfield   = &(field[numls]);
  numff       = genprob.numff;
  fluidfield  = &(field[numff]);
  /* plausibility check */
  dsassert(lsetfield->fieldtyp==levelset,"FIELD 0 has to be levelset\n");
  dsassert(fluidfield->fieldtyp==fluid,"FIELD 1 has to be fluid\n");

  /* I N I T I A L I Z A T I O N */
  eflag=1;
  /* set */
  lsdyn= alldyn[numls].lsdyn;
  lsdata = lsdyn->lsdata;
  fdyn = alldyn[numff].fdyn;
  /* init all applied time curves */
  for (actcurve = 0;actcurve<numcurve;actcurve++)
    dyn_init_curve(actcurve,lsdyn->nstep,lsdyn->dt,lsdyn->maxtime);
  /* initialize fluid */
  ls_fluid(eflag);
  /* initialise levelset */
  ls_levelset(eflag);
  /* initialize level set module */
  lsflag = ls_initphase;
  ls_main(lsflag);
  
  /* T I M E L O O P */
 timeloop:
  eflag=2;
  lsdyn->step++;
  lsdyn->time += lsdyn->dt; 
  fdyn->step = lsdyn->step;
  fdyn->acttime = lsdyn->time;
  
  /* solve fluid field */
  ls_fluid(eflag);
  /* solve levelset field */
  if (lsdyn->step>lsdyn->nfstep)
  {
    ls_levelset(eflag);
  }

  /* F I N A L I Z E */
  eflag=3;
  /* finalize fluid field */
  ls_fluid(eflag);
  /* finalize levelset field */
  if (lsdyn->step>lsdyn->nfstep)
  {
    ls_levelset(eflag);
  }
  /* perform reinitialization */
  if (lsdyn->step>lsdyn->nfstep)
  {
    if (lsdyn->step%lsdyn->nfreinit==0)
    {
      if (lsdyn->lsdata->reinitialization==1)
      {
        /* copy solution from sol_increment[1][j] to sol_increment[2][j] */
        solserv_sol_copy(lsetfield,0,1,1,1,2);
        /* turn reinitialization flag on */
        if (lsdyn->lsdata->reinitflag==1) dserror("reinitflag is on!");
        lsdyn->lsdata->reinitflag = 1;
        printf("\n**WARNING** PERFORMING RE-INITIALIZATION!");
        for (i=0; i<lsdyn->lsdata->numreinit; i++)
        {
          /* solve levelset field */
          eflag = 2;
          ls_levelset(eflag);
          /* finalize levelset field */      
          eflag = 3;
          ls_levelset(eflag);
        }
        /* turn reinitialization flag off */
        lsdyn->lsdata->reinitflag = 0;
      }
    }
  }
  /* update level set module */
  if (lsdyn->step>lsdyn->nfstep)
  {
    lsflag = ls_updtphase;
    ls_main(lsflag);
  }
  /* finalising this time step */
  if (lsdyn->step < lsdyn->nstep && lsdyn->time <= lsdyn->maxtime)
    goto timeloop;
  
/* C L E A N I N G   U P   P H A S E */
  eflag=99;
  ls_fluid(eflag);
  ls_levelset(eflag);
  /* finalize level set module */
  lsflag = ls_finaphase;
  ls_main(lsflag);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  return;
} /* end of ls_dyn */
#endif
