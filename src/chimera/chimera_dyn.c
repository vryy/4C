#ifdef D_CHIMERA
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "chimera.h"



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
extern struct _PARTITION        *partition;
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
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS          ioflags;
/*----------------------------------------------------------------------*
  | global variable *solv, vector of lenght numfld of structures SOLVAR  |
  | defined in solver_control.c                                          |
  |                                                                      |
  |                                                       m.gee 11/00    |
  *----------------------------------------------------------------------*/
extern struct _SOLVAR           *solv;
/*-----------------------------------------------------------------------*
  | enum _CALC_ACTION                                      m.gee 1/02    |
  | command passed from control routine to the element level             |
  | to tell element routines what to do                                  |
  | defined globally in global_calelm.c                                  |
  *----------------------------------------------------------------------*/
extern enum   _CALC_ACTION       calc_action[MAXFIELD];
extern struct _CHIMERA_DATA     *chm_data;





void chimera_dyn()
{
  INT                i,rr;
  FIELD             *fluidfield;
  FLUID_DYNAMIC     *fdyn;

  INTRA             *actintra;          /* pointer to active intra-comm. */
  PARTITION         *actpart;           /* pointer to active partition   */
  CALC_ACTION       *action;

  INT                numdis;
  INT                Konvergenz_Gebietsiteration;
  INT                nmax_Gebietsiteration;
  INT                steady = 0;
  INT                actpos = 0;

#ifdef DEBUG
  dstrc_enter("chimera_dyn");
#endif
/*----------------------------------------------------------------------*/

  /* plausibility check */
  dsassert(genprob.numfld==1,"**ERROR** numfld!=1 for chimera problem\n");
  dsassert(genprob.numff==0,"**ERROR** genprob.numff!=0 for chimera problem\n");

  /* access to the fluid field */
  fluidfield = &(field[genprob.numff]);
  fdyn = alldyn[genprob.numff].fdyn;
  action  = &(calc_action[genprob.numff]);

  /* if we are not parallel, we have to allocate an alibi
     intra-communicator structure */
  actpart     = &(partition[genprob.numff]);

#ifdef PARALLEL
  actintra    = &(par.intra[genprob.numff]);
#else
  actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
  if (!actintra) dserror("Allocation of INTRA failed");
  actintra->intra_fieldtyp = fluid;
  actintra->intra_rank     = 0;
  actintra->intra_nprocs   = 1;
#endif

  /* init all applied time curves */
  for (i=0; i<numcurve; i++)
    dyn_init_curve(i,fdyn->nstep,fdyn->dt,fdyn->maxtime);

  if (chm_data->interact_dis_on_off==1)
  {
/****************************CHIMERA REGION START************************/
/****************************CHIMERA REGION START************************/
/****************************CHIMERA REGION START************************/

    /* initialize quadtree structure */
    if (chm_data->search_action==1)
    {
      for (i=0; i<fluidfield->ndis; i++)
      {
      chimera_search_quadtree_init(&(fluidfield->dis[i]),i,chm_data->nmax_pro_Blatt_in_quadtree);
      }
    }

    /* perform automatic hole cutting */
    /*
      NOTE =>
      at the moment there is no need to perform the initialiZation in
      every timestep <---> moving body problems have not been yet considered
      (vorerst nicht in jedem Zeitschritt, da ohne Bewegung des Gitters)
    */
    if (chm_data->chimera_define_back_bound==define_background_boundary_automatic)
    {
      chimera_automatic_hole_cutting(chm_data->n_Hintergrunddiskretisierung);
    }

/****************************CHIMERA REGION END**************************/
/****************************CHIMERA REGION END**************************/
/****************************CHIMERA REGION END**************************/
  }

  /* write mesh data into gid file */
  if (ioflags.fluid_sol_gid==1 && par.myrank==0)
  {
    out_gid_msh_trial();
  }

  /* initialize the fluid problem on overlapping domains */
  /*
    NOTE =>
    has to be called after chimera_automatic_hole_cutting to get correct
    initialisation for sparse arrays
    mu"s nach chimera_automatic_hole_cutting kommen, sonst werden
    sparse-arrays falsch initialisiert
  */
  chimera_initialize();
#if 0
  if (ioflags.chimera_to_matlab==1)
  {
    chimera_to_matlab();
  }
#endif

  /* print out initial data to .out */
  out_sol(fluidfield,actpart,actintra,fdyn->step,actpos);

  if (ioflags.fluid_sol_gid==1 && par.myrank==0)
  {
    out_gid_sol_trial("velocity",fluidfield,actintra,fdyn->step,actpos,fdyn->acttime);
    out_gid_sol_trial("pressure",fluidfield,actintra,fdyn->step,actpos,fdyn->acttime);
  }

  /************************************************************************/
  /*                         START TIME LOOP                              */
  /************************************************************************/
 timeloop:
  fdyn->step++;
  fdyn->acttime += fdyn->dt;

  if (par.myrank==0)
  {
    printf("\n\n\n");
    printf("***************************************************************\n\n");
    printf("TIMESTEP %d   (INTEGRATION SCHEME => ONE-STEP THETA)\n",fdyn->step);
    printf("\n\n***************************************************************\n\n");
  }

  if (chm_data->interact_dis_on_off==1)
  {
/****************************CHIMERA REGION START************************/
/****************************CHIMERA REGION START************************/
/****************************CHIMERA REGION START************************/

    /* Sichern der Information aus dem Zeitschritt f"ur Konvergenz"uberpr"ufung */
    /* copy solution from sol_increment[3][j] to sol_increment[5][j] */
    for (rr=0; rr<fluidfield->ndis; rr++)
    {
      solserv_sol_copy(fluidfield,rr,1,1,3,5);
    }

/****************************CHIMERA REGION END**************************/
/****************************CHIMERA REGION END**************************/
/****************************CHIMERA REGION END**************************/
  }

  if (chm_data->interact_dis_on_off==1)
  {
/****************************CHIMERA REGION START************************/
/****************************CHIMERA REGION START************************/
/****************************CHIMERA REGION START************************/

    /**********************************************************************/
    /*                     START SUBDOMAIN ITERATION                      */
    /**********************************************************************/
    chm_data->Konvergenz_Gebietsiteration = 0;        /*  Schalter  */
    rr = 0;

    while (rr<chm_data->nmax_Gebietsiteration && chm_data->Konvergenz_Gebietsiteration==0)
    {
      if (par.myrank==0)
      {
        printf("\n\n*****************************************************************************\n");
        printf("SUBDOMAIN ITERATION %3d",rr+1);
      }

      chm_data->rel_err[0]=0;
      chm_data->rel_err[1]=0;

      chm_data->abs_err[0]=0;
      chm_data->abs_err[1]=0;

      /* solve fluid problem on each discretization */
      for (i=0; i<fluidfield->ndis; i++) /* => i == numdis */
      {
        if (par.myrank==0)
        {
          printf("\n\n**WARNING** SOLVING SUBDOMAIN [%d]\n",i+1);
        }
        /* solve time step till convergence */
        chimera_solve(i);
      }
      /* erh"ohen des Iterationsz"ahlers*/
      /* increment iteration counter */
      rr++;

      /* convergence of subdomain iteration*/
      if (par.myrank==0)
      {
        printf("=> CONVERGENCE CHECK FOR SUBDOMAIN ITERATION\n\n");
        printf("                  SUBDOMAIN 1             SUBDOMAIN 2\n");
        printf("   REL_ERR:      %10.7E           %10.7E\n",chm_data->rel_err[0],chm_data->rel_err[1]);
        printf("   ABS_ERR:      %10.7E           %10.7E\n",chm_data->abs_err[0],chm_data->abs_err[1]);
        printf("\n*****************************************************************************\n");
      }
      /* Abbruch bei unterschreiten der Toleranz. Mindestens zwei
         Iterationsschritte */
      /* break loop if current change in boundary value is lower then a
         given tolerance 2 Iterations minimum */
      if (rr>1)
      {
        if ((chm_data->rel_err[0]<chm_data->rel_TOL &&
             chm_data->rel_err[1]<chm_data->rel_TOL)||
            (chm_data->abs_err[0]<chm_data->abs_TOL &&
             chm_data->abs_err[1]<chm_data->abs_TOL))
        {
          chm_data->Konvergenz_Gebietsiteration=1;
        }
      }
    }
    /************************************************************************/
    /*                        END SUBDOMAIN ITERATION                       */
    /************************************************************************/

    /* finalize iteration step */
    for (i=0; i<fluidfield->ndis; i++)
    {
      chimera_finalize(i);
      /* copy solution from sol_increment to sol */
      fluid_sol_copy(fluidfield,i,1,0,1,0,fdyn->numdf);
    }

    /* perform steady check */
    steady = fluid_steadycheck_chimera(fluidfield,solv[0].sol[0].numeq_total);
    /* write solution to .flavia.res */
    if (ioflags.fluid_sol_gid==1 && par.myrank==0)
    {
      if (fdyn->step%chm_data->write_out_every==0)
      {
        out_gid_sol_trial("velocity",fluidfield,actintra,fdyn->step,0,fdyn->acttime);
        out_gid_sol_trial("pressure",fluidfield,actintra,fdyn->step,0,fdyn->acttime);
      }
    }

    /* Ausgabe des aktuellen Ergebnisses */
    /* mit pager */
    if (fdyn->step%chm_data->write_out_every==0)
    {
      if (ioflags.chimera_to_pager==1 && par.myrank==0)
      {
        for (numdis=0; numdis<fluidfield->ndis; numdis++)
        {
          chimera_output(numdis);
        }
      }
    }

    /* write solution to matlab */
#if 0
    if (ioflags.chimera_to_matlab==1)
    {
      chimera_write_soln_to_matlab();
    }
#endif

/****************************CHIMERA REGION END**************************/
/****************************CHIMERA REGION END**************************/
/****************************CHIMERA REGION END**************************/

  }
  else if (chm_data->interact_dis_on_off==0)
  {
    /* solve fluid problem on each discretization */
    for (i=0; i<fluidfield->ndis; i++)               /* => i == numdis  */
    {
      /* solve time step till convergence */
      chimera_solve(i);
      /* finalize time step */
      chimera_finalize(i);
    }

    /* write solution to matlab */
#if 0
    if (ioflags.chimera_to_matlab==1)
    {
      chimera_write_soln_to_matlab();
    }
#endif
  }

  /* check time step */
  if (fdyn->step < fdyn->nstep && fdyn->acttime <= fdyn->maxtime)
    goto timeloop;
  /************************************************************************/
  /*                           END TIME LOOP                              */
  /************************************************************************/

  /* Ausgabe des Endergebnisses */
  /* mit pager */
  if (ioflags.chimera_to_pager==1 && par.myrank==0)
  {
    for (numdis=0; numdis<fluidfield->ndis; numdis++){
      chimera_output(numdis);
      printf("ne %5d\nnn %5d\n",field[0].dis[numdis].numele,field[0].dis[numdis].numnp);
      printf("boundary on\ncolor on\nncolor %d\n",numdis+1);
      printf("nen %d\n",field[0].dis[numdis].element[0].numnp);
      printf("mien ../output/chimera_ien%d\n",numdis);
      printf("mxyz ../output/chimera_xy%d\n",numdis);
      printf("data ../output/chimera_data%d\n",numdis);
      if(numdis == 0)
      {
        printf("line_width 0.2\n");
        printf("vect_size 0.025\n");
        printf("vect_scale 0.05\n");
        printf("little-endian on\n");
        printf("ndf 3\n");
        printf("nrec 1\n");
        printf("idf 3\n");
        printf("xmax 1.5\n");
        printf("ymax 1.5\n");
        printf("xmin -0.5\n");
        printf("ymin -0.5\n");
      }
      printf("vector\n");
      printf("plot\n");
    }
  }

  /* clean up */
  chimera_clean();

  if (chm_data->interact_dis_on_off==1)
  {
/****************************CHIMERA REGION START************************/
/****************************CHIMERA REGION START************************/
/****************************CHIMERA REGION START************************/

    if (chm_data->search_action==1)
    {
      for(rr=0;rr<2;rr++)
      {
        chimera_search_quadtree_free(rr);
      }
    }
    free(chm_data->rel_err);
    free(chm_data->abs_err);
    free(chm_data);

/****************************CHIMERA REGION END**************************/
/****************************CHIMERA REGION END**************************/
/****************************CHIMERA REGION END**************************/
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of chimera_dyn */
#endif /* D_CHIMERA */
