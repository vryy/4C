/*!----------------------------------------------------------------------
\file
\brief ls_dyn.c

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_LS
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../xfem/xfem_prototypes.h"
#include "ls_prototypes.h"
/*!
\addtogroup LEVELSET
*//*! @{ (documentation module open)*/



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA           *alldyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB    genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD     *field;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR        par;
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern INT                numcurve;
extern struct _CURVE     *curve;



/*!----------------------------------------------------------------------
\brief control subroutine for level set problem

<pre>                                                            irhan 05/04
in this subroutine a time history analysis is performed for two phase fluid
problems. the method used is coupled level set / extended finite element
method.
</pre>

*----------------------------------------------------------------------*/
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
  INT                pcnt = 0;

  LS_GEN_DATA       *lsdata;
  LSFLAG             lsflag;

#ifdef DEBUG
  dstrc_enter("ls_dyn");
#endif
/*----------------------------------------------------------------------*/

  switch (genprob.probtyp)
  {
      case prb_twophase:
        /*
         * convention used at the moment =>
         * TWO PHASE FLUID FLOW -> FIELD 0: levelset
         *                         FIELD 1: fluid
         */

        /* I N I T I A L I Z A T I O N */
        /* access to field objects */
        numfld=genprob.numfld;
        dsassert(numfld==2,"Two fields needed for Two Phase Fluid Flow Problem!\n");
        numls       = genprob.numls;
        lsetfield   = &(field[numls]);
        numff       = genprob.numff;
        fluidfield  = &(field[numff]);
        /* plausibility check */
        dsassert(lsetfield->fieldtyp==levelset,"FIELD 0 has to be levelset\n");
        dsassert(fluidfield->fieldtyp==fluid,"FIELD 1 has to be fluid\n");
        /* set */
        lsdyn= alldyn[numls].lsdyn;
        lsdata = lsdyn->lsdata;
        fdyn = alldyn[numff].fdyn;
        /* set execution flag */
        eflag=1;
        /* check */
        if (lsdata->setvel!=2)
          dserror("lsdata->setvel!=2 for Two Phase Fluid Flow Problem");
        /* init all applied time curves */
        for (actcurve = 0;actcurve<numcurve;actcurve++)
          dyn_init_curve(actcurve,lsdyn->nstep,lsdyn->dt,lsdyn->maxtime);
        /* initialize fluid problem */
        ls_fluid(eflag);
        /* initialize levelset problem */
        ls_levelset(eflag);
        /* initialize levelset module */
        lsflag = ls_initphase;
        ls_main(lsflag);
        /* T I M E L O O P */
  timeloop1:
        lsdyn->step++;
        lsdyn->time += lsdyn->dt;
        fdyn->step = lsdyn->step;
        fdyn->acttime = lsdyn->time;
        /* set execution flag */
        eflag=2;
        /* solve fluid field */
        ls_fluid(eflag);
        /* solve levelset field */
        if (lsdyn->step>lsdyn->nfstep)
          ls_levelset(eflag);
        /* F I N A L I Z E */
        /* set execution flag */
        eflag=3;
        /* finalize fluid field */
        ls_fluid(eflag);
        /* finalize levelset field */
        if (lsdyn->step>lsdyn->nfstep)
          ls_levelset(eflag);
        /* START REINITIALIZATION */
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
              /* print to screen */
              printf("\n**WARNING** PERFORMING RE-INITIALIZATION!");
              if (lsdyn->lsdata->anchor == 1)
              {
                printf("\n**WARNING** INTERFACE IS ANCHORED!");
              }
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
        } /* END REINITIALIZATION */
        /* update level set module */
        if (lsdyn->step>lsdyn->nfstep)
        {
          /* increment the print counter (used to print intermediate results into files) */
          pcnt++;
          if (pcnt==lsdyn->nfreprn)
          {
            if (lsdyn->lsdata->print_on_off!=0) dserror("lsdyn->lsdata->print_on_off!=0");
            /* turn on print flag */
            lsdyn->lsdata->print_on_off = 1;
            pcnt = 0;
          }
          lsflag = ls_updtphase;
          ls_main(lsflag);
        }
        /* finalizing this time step */
        if (lsdyn->step < lsdyn->nstep && lsdyn->time <= lsdyn->maxtime)
          goto timeloop1;
        /* C L E A N I N G   U P   P H A S E */
        /* set execution flag */
        eflag=99;
        ls_fluid(eflag);
        ls_levelset(eflag);
        /* finalize level set module */
        lsflag = ls_finaphase;
        ls_main(lsflag);
        break;
      case prb_levelset:
        /*
         * convention used at the moment =>
         * PURE LEVELSET PROBLEM -> FIELD 0: levelset
         */

        /* I N I T I A L I Z A T I O N */
        /* access to field object */
        numfld=genprob.numfld;
        dsassert(numfld==1,"One field needed for Pure Levelset Problem!\n");
        numls       = genprob.numls;
        lsetfield   = &(field[numls]);
        /* plausibility check */
        dsassert(lsetfield->fieldtyp==levelset,"FIELD 0 has to be levelset\n");
        /* set */
        lsdyn= alldyn[numls].lsdyn;
        lsdata = lsdyn->lsdata;
        /* set execution flag */
        eflag=1;
        /* check */
        if (lsdata->setvel!=1)
          dserror("lsdata->setvel!=1 for Pure Levelset Problem");
        /* init all applied time curves */
        for (actcurve = 0;actcurve<numcurve;actcurve++)
          dyn_init_curve(actcurve,lsdyn->nstep,lsdyn->dt,lsdyn->maxtime);
        /* initialize levelset problem */
        ls_levelset(eflag);
        /* initialize levelset module */
        lsflag = ls_initphase;
        ls_main(lsflag);
        /* T I M E L O O P */
  timeloop2:
        lsdyn->step++;
        lsdyn->time += lsdyn->dt;
        /* set execution flag */
        eflag=2;
        /* solve levelset field */
        ls_levelset(eflag);
        /* F I N A L I Z E */
        /* set execution flag */
        eflag=3;
        /* finalize levelset field */
        ls_levelset(eflag);
        /* START REINITIALIZATION */
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
              /* print to screen */
              printf("\n**WARNING** PERFORMING RE-INITIALIZATION!");
              if (lsdyn->lsdata->anchor == 1)
              {
                printf("\n**WARNING** INTERFACE IS ANCHORED!");
              }
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
        } /* END REINITIALIZATION */
        /* increment the print counter (used to print intermediate results into files) */
        pcnt++;
        if (pcnt==lsdyn->nfreprn)
        {
          if (lsdyn->lsdata->print_on_off!=0) dserror("lsdyn->lsdata->print_on_off!=0");
          /* turn on print flag */
          lsdyn->lsdata->print_on_off = 1;
          pcnt = 0;
        }
        /* update level set module */
        lsflag = ls_updtphase;
        ls_main(lsflag);
        /* finalizing this time step */
        if (lsdyn->step < lsdyn->nstep && lsdyn->time <= lsdyn->maxtime)
          goto timeloop2;
        /* C L E A N I N G   U P   P H A S E */
        /* set execution flag */
        eflag=99;
        ls_levelset(eflag);
        /* finalize level set module */
        lsflag = ls_finaphase;
        ls_main(lsflag);
        break;
      default:
        dserror("in ls_dyn() unknown Problemtyp requested");
        break;
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of ls_dyn */
/*! @} (documentation module close)*/
#endif
