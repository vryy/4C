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
extern INT                 numcurve;
extern struct _CURVE      *curve;
extern struct _LS_DYNAMIC *lsdyn;
extern struct _XFEM_DATA   xfem_data;


struct _FIELD      *fluidfield;
struct _FIELD      *lsetfield;


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
  FLUID_DYNAMIC     *fdyn;
  INT                numfld;
  INT                numff;
  INT                numls;
  INT                actcurve;
  INT                datastep = 0;     /* counter for output print */
  INT                restartstep = 0;  /* counter for restart data print */


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
         * TWO PHASE FLUID FLOW -> FIELD 0: fluid
         *                         FIELD 1: levelset
         */

        /* I N I T I A L I Z A T I O N */
        /* access to field objects */
        numfld = genprob.numfld;
        dsassert(numfld==2,"Two fields needed for Two Phase Fluid Flow Problem!\n");
        numls = genprob.numls;
        lsetfield = &(field[numls]);
        numff = genprob.numff;
        fluidfield = &(field[numff]);
        /* plausibility check */
        dsassert(lsetfield->fieldtyp==levelset,"FIELD 0 has to be levelset\n");
        dsassert(fluidfield->fieldtyp==fluid,"FIELD 1 has to be fluid\n");
        /* set */
        lsdata = lsdyn->lsdata;
        fdyn = alldyn[numff].fdyn;
        /* set execution flag */
        eflag = 1;
        /* check */
        if (lsdata->setvel != 2)
          dserror("lsdata->setvel!=2 for Two Phase Fluid Flow Problem");
        /* init all applied time curves */
        for (actcurve = 0;actcurve<numcurve;actcurve++)
          dyn_init_curve(actcurve,lsdyn->nstep,lsdyn->dt,lsdyn->maxtime);

        if (xfem_data.xfem_optimize == 1)
        {
          /* init degrees of freedom for levelset field */
          init_dof_discretization(&(lsetfield->dis[0]));
        }

        /* initialize levelset problem */
        ls_levelset(eflag);

        /* initialize levelset module */
        lsflag = ls_initphase;
        ls_main(lsflag);

        if (xfem_data.xfem_optimize == 1)
        {
          /* init degrees of freedom for fluid field */
          init_dof_discretization_xfem(&(fluidfield->dis[0]));
        }

        /* initialize fluid problem */
        ls_fluid(eflag);

        if (xfem_data.xfem_optimize == 1)
        {
          /* assign degrees of freedom to level set field */
          assign_dof_discretization(&(lsetfield->dis[0]));
          /* assign degrees of freedom to fluid field */
          assign_dof_discretization_xfem(&(fluidfield->dis[0]));
          /* mask matrices */
          mask_global_matrices();
        }
        eflag = 4;
        /* initialize fluid solver */
        ls_fluid(eflag);
        /* initialize levelset solver */
        ls_levelset(eflag);

        /* T I M E L O O P */
  timeloop1:
        lsdyn->step++;
        lsdyn->acttime += lsdyn->dt;
        fdyn->step = lsdyn->step;
        fdyn->acttime = lsdyn->acttime;

        /* set execution flag */
        eflag = 2;

        /* solve fluid field */
        ls_fluid(eflag);

        /* solve levelset field */
        if (lsdyn->step > lsdyn->nfstep)
        {
          /* set flag for element calculation */
          lsdyn->lsdata->ls2_calc_flag = ls2_advection;

          ls_levelset(eflag);
        }

        /* F I N A L I Z E */

        /* set execution flag */
        eflag = 3;

        /* finalize fluid field */
        ls_fluid(eflag);

        /* finalize levelset field */
        if (lsdyn->step > lsdyn->nfstep)
        {
          ls_levelset(eflag);
        }

        /* update level set module */
        if (lsdyn->step > lsdyn->nfstep)
        {
          lsflag = ls_updtphase;
          ls_main(lsflag);
        }

        /* START REINITIALIZATION */
        if (lsdyn->step > lsdyn->nfstep)
        {
          if (lsdyn->step%lsdyn->perf_reinit_evry == 0)
          {
            if (lsdyn->lsdata->reinitialization == 1)
            {
              /* set flag for element calculation */
              lsdyn->lsdata->ls2_calc_flag = ls2_reinitialization;

              /* copy solution from sol_increment[1][j] to sol_increment[2][j] */
              /* this operation is necessary for the computation of sign function */
              solserv_sol_copy(lsetfield,0,node_array_sol_increment,node_array_sol_increment,1,2);

              /* print to screen */
              printf("\n**WARNING** PERFORMING RE-INITIALIZATION!");
              if (lsdyn->lsdata->anchor == 1)
              {
                printf("\n**WARNING** INTERFACE IS ANCHORED!");
              }
              for (i=0; i<lsdyn->lsdata->numreinit; i++)
              {
                /* set execution flag */
                eflag = 2;

                /* solve levelset field */
                ls_levelset(eflag);

                /* set execution flag */
                eflag = 3;

                /* finalize levelset field */
                ls_levelset(eflag);
              }
            }
          }
        }
        /* END REINITIALIZATION */






/****************** TEST VERSION *************************/
/****************** TEST VERSION *************************/
/****************** TEST VERSION *************************/
/****************** TEST VERSION *************************/
/****************** TEST VERSION *************************/
/****************** TEST VERSION *************************/




        if (lsdyn->step > lsdyn->nfstep)
        {

          /* START COMPUTING GRADIENT FIELD */

          /* set flag for element calculation */
          lsdyn->lsdata->ls2_calc_flag = ls2_gradient;

          printf("\n**WARNING** COMPUTING GRADIENT FIELD!");

          /* COMPUTE FIRST COMPONENT OF THE GRADIENT FIELD */

          /* set spatial direction ( 1 or 2 ) */
          lsdyn->lsdata->gradient_direction = 1;

          /* set execution flag */
          eflag = 2;

          /* solve levelset field */
          ls_levelset(eflag);



          /* COMPUTE SECOND COMPONENT OF THE GRADIENT FIELD */

          /* set spatial direction ( 1 or 2 ) */
          lsdyn->lsdata->gradient_direction = 2;

          /* set execution flag */
          eflag = 2;

          /* solve levelset field */
          ls_levelset(eflag);

          /* END COMPUTING GRADIENT FIELD */
        }




/****************** TEST VERSION *************************/
/****************** TEST VERSION *************************/
/****************** TEST VERSION *************************/
/****************** TEST VERSION *************************/
/****************** TEST VERSION *************************/
/****************** TEST VERSION *************************/









        /* write data to files */
        if (lsdyn->step > lsdyn->nfstep)
        {
          /* increment the print counters (used to print intermediate results into files) */
          datastep++;
          restartstep++;
          if (datastep == lsdyn->data_write_evry)
          {
            if (lsdyn->lsdata->print_on_off != 0) dserror("lsdyn->lsdata->print_on_off!=0");
            /* turn on print flag */
            lsdyn->lsdata->print_on_off = 1;
            datastep = 0;
          }
          lsflag = ls_writphase;
          ls_main(lsflag);

          /* write restart to pss file */
          if (restartstep == lsdyn->res_write_evry)
          {
            restartstep=0;
            eflag = 5;
            ls_levelset(eflag);
          }
        }

        if (xfem_data.xfem_optimize == 1)
        {
          /* assign degrees of freedom to fluid field */
          assign_dof_discretization_xfem(&(fluidfield->dis[0]));
          /* mask matrices */
          remask_global_matrices_xfem();
        }
        eflag = 5;
        /* reinitialize fluid solver */
        ls_fluid(eflag);

        /* finalizing this time step */
        if (lsdyn->step < lsdyn->nstep && lsdyn->acttime <= lsdyn->maxtime)
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
	if (lsdyn->lsdata->probdescr == 3) /* advection equation is solved */
	  {
	    /*
	     * convention used at the moment =>
	     * PURE LEVELSET PROBLEM -> FIELD 0: levelset
	     */
	    /* I N I T I A L I Z A T I O N */
	    /* access to field objects */
	    numfld = genprob.numfld;
	    dsassert(numfld==1,"Two fields needed for Pure Level Set Problem!\n");
	    numls       = genprob.numls;
	    lsetfield   = &(field[numls]);
	    /* plausibility check */
	    dsassert(lsetfield->fieldtyp==levelset,"FIELD 0 has to be levelset\n");
	    /* set */
	    lsdata = lsdyn->lsdata;
	    /* set execution flag */
	    eflag=1;
	    /* check */
	    if (lsdata->setvel != 1)
	      dserror("lsdata->setvel!=1 for Pure Level Set Problem");
	    /* init all applied time curves */
	    for (actcurve = 0;actcurve<numcurve;actcurve++)
	      dyn_init_curve(actcurve,lsdyn->nstep,lsdyn->dt,lsdyn->maxtime);

            if (xfem_data.xfem_optimize == 1)
            {
              /* init degrees of freedom for levelset field */
              init_dof_discretization(&(lsetfield->dis[0]));
            }

            /* initialize levelset problem */
	    ls_levelset(eflag);

            if (xfem_data.xfem_optimize == 1)
            {
              /* assign degrees of freedom to level set field */
              assign_dof_discretization(&(lsetfield->dis[0]));
              /* mask matrices */
              mask_global_matrices();
            }
            eflag = 4;
            /* initialize levelset solver */
            ls_levelset(eflag);


            printf("i am here!\n");

	    /* initialize levelset module */
	    lsflag = ls_initphase;
	    ls_main(lsflag);
	    /* T I M E L O O P */
	  timeloop2:
	    lsdyn->step++;
	    lsdyn->acttime += lsdyn->dt;
	    /* set execution flag */
	    eflag=2;
            /* set flag for element calculation */
            lsdyn->lsdata->ls2_calc_flag = ls2_advection;
	    /* solve levelset field */
	    ls_levelset(eflag);
	    /* F I N A L I Z E */
	    /* set execution flag */
	    eflag=3;
	    /* finalize levelset field */
	    ls_levelset(eflag);
	    /* update level set module */
	    lsflag = ls_updtphase;
	    ls_main(lsflag);
	    /* START REINITIALIZATION */
	    if (lsdyn->step%lsdyn->perf_reinit_evry == 0)
            {
              if (lsdyn->lsdata->reinitialization == 1)
              {
                /* set flag for element calculation */
                lsdyn->lsdata->ls2_calc_flag = ls2_reinitialization;
                /* copy solution from sol_increment[1][j] to sol_increment[2][j] */
                solserv_sol_copy(lsetfield,0,1,1,1,2);
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
              }
            } /* END REINITIALIZATION */
	    /* increment the print counter (used to print intermediate results into files) */
	    datastep++;
	    restartstep++;
	    if (datastep == lsdyn->data_write_evry)
            {
              if (lsdyn->lsdata->print_on_off != 0) dserror("lsdyn->lsdata->print_on_off!=0");
              /* turn on print flag */
              lsdyn->lsdata->print_on_off = 1;
              datastep = 0;
            }
            /* write data to files */
	    lsflag = ls_writphase;
	    ls_main(lsflag);
	    /* write restart to pss file */
	    if (restartstep == lsdyn->res_write_evry)
	      {
		restartstep=0;
		eflag = 5;
		ls_levelset(eflag);
	      }
	    /* finalizing this time step */
	    if (lsdyn->step < lsdyn->nstep && lsdyn->acttime <= lsdyn->maxtime)
	      goto timeloop2;
	    /* C L E A N I N G   U P   P H A S E */
	    /* set execution flag */
	    eflag=99;
	    ls_levelset(eflag);
	    /* finalize level set module */
	    lsflag = ls_finaphase;
	    ls_main(lsflag);
	  }
	else if (lsdyn->lsdata->probdescr == 4) /* reinitialization is performed */
	  {
	    /*
	     * convention used at the moment =>
	     * PURE LEVELSET PROBLEM -> FIELD 0: levelset
	     */

	    /* I N I T I A L I Z A T I O N */
	    /* access to field objects */
	    numfld = genprob.numfld;
	    dsassert(numfld==1,"Two fields needed for Pure Level Set Problem!\n");
	    numls       = genprob.numls;
	    lsetfield   = &(field[numls]);
	    /* plausibility check */
	    dsassert(lsetfield->fieldtyp==levelset,"FIELD 0 has to be levelset\n");
	    /* set */
	    lsdata = lsdyn->lsdata;
	    /* set execution flag */
	    eflag=1;
	    /* check */
	    if (lsdata->setvel != 1)
	      dserror("lsdata->setvel!=1 for Pure Level Set Problem");
	    /* init all applied time curves */
	    for (actcurve = 0;actcurve<numcurve;actcurve++)
	      dyn_init_curve(actcurve,lsdyn->nstep,lsdyn->dt,lsdyn->maxtime);
	    /* initialize levelset problem */
	    ls_levelset(eflag);
            eflag = 4;
            /* initialize levelset solver */
            ls_levelset(eflag);
	    /* initialize levelset module */
	    lsflag = ls_initphase;
	    ls_main(lsflag);
	    /* copy solution from sol_increment[1][j] to sol_increment[2][j] */
	    solserv_sol_copy(lsetfield,0,1,1,1,2);
	    /* set */
	    lsdyn->dt = lsdyn->lsdata->rdt;
            lsdyn->nstep = lsdyn->lsdata->numreinit;

	    /* T I M E L O O P */
  timeloop3:
	    lsdyn->step++;
	    lsdyn->acttime += lsdyn->dt;

            /* set flag for element calculation */
            lsdyn->lsdata->ls2_calc_flag = ls2_reinitialization;
	    /* print to screen */
	    printf("\n**WARNING** PERFORMING RE-INITIALIZATION!");
	    if (lsdyn->lsdata->anchor == 1)
	      {
		printf("\n**WARNING** INTERFACE IS ANCHORED!");
	      }
	    /* solve levelset field */
	    eflag = 2;
	    ls_levelset(eflag);
	    /* finalize levelset field */
	    eflag = 3;
	    ls_levelset(eflag);

	    /* increment the print counter (used to print intermediate results into files) */
	    datastep++;
	    restartstep++;
	    if (datastep == lsdyn->data_write_evry)
            {
		if (lsdyn->lsdata->print_on_off != 0) dserror("lsdyn->lsdata->print_on_off!=0");
		/* turn on print flag */
		lsdyn->lsdata->print_on_off = 1;
		datastep = 0;
            }
	    lsflag = ls_writphase;
	    ls_main(lsflag);

	    /* write restart to pss file */
	    if (restartstep == lsdyn->res_write_evry)
            {
              restartstep=0;
              eflag = 5;
              ls_levelset(eflag);
            }

            /* finalizing this time step */
	    if (lsdyn->step < lsdyn->nstep && lsdyn->acttime <= lsdyn->maxtime)
	      goto timeloop3;
	    /* C L E A N I N G   U P   P H A S E */
	    /* set execution flag */
	    eflag=99;
	    ls_levelset(eflag);
	    /* finalize level set module */
	    lsflag = ls_finaphase;
	    ls_main(lsflag);
	  }
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
