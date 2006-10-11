/*!----------------------------------------------------------------------
\file
\brief input of control information and general problem data

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
struct _FILES           allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA       *alldyn;                                               |
 *----------------------------------------------------------------------*/
extern ALLDYNA        *alldyn;
/*----------------------------------------------------------------------*
 |                                                          al 08/02    |
 | pointer to allocate eigensolution variables                          |
 | dedfined in global_control.c                                         |
 | struct _ALLEIG       *alleig;                                        |
 *----------------------------------------------------------------------*/
extern ALLEIG              *alleig;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _STATIC_VAR  *statvar;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 | struct _FIELD         *field;                                        |
 *----------------------------------------------------------------------*/
extern struct _FIELD       *field;

/*----------------------------------------------------------------------*
 | input of control information                           m.gee 8/00    |
 *----------------------------------------------------------------------*/


#ifdef WALLCONTACT
extern struct _WALL_CONTACT contact;
#endif

#ifdef D_SSI
void inpctr_dyn_ssi(SSI_DYNAMIC *ssidyn);
#endif

void inpctr()
{
#ifdef DEBUG
dstrc_enter("inpctr");
#endif

/*---------------------------------------------- input of problem title */
   inpctrhed();
/*--------------------------------------- input of general problem data */
   inpctrprob();
/*---------------------------------------- input of general solver data */
   /* for FSI */
   if (genprob.probtyp == prb_fsi)
   {
      if (genprob.numfld!=3) dserror("numfld != 3 for FSI");

      solv = (SOLVAR*)CCACALLOC(genprob.numfld,sizeof(SOLVAR));

      solv[genprob.numsf].fieldtyp = structure;
      inpctrsol(&(solv[genprob.numsf]));

      solv[genprob.numff].fieldtyp = fluid;
      inpctrsol(&(solv[genprob.numff]));

      solv[genprob.numaf].fieldtyp = ale;
      inpctrsol(&(solv[genprob.numaf]));
   }
   /* for SSI */
   if (genprob.probtyp == prb_ssi)
   {
      if (genprob.numfld!=2) dserror("numfld != 2 for SSI");

      solv = (SOLVAR*)CCACALLOC(genprob.numfld,sizeof(SOLVAR));

      solv[0].fieldtyp = structure;
      inpctrsol(&(solv[0]));

      solv[1].fieldtyp = structure;
      inpctrsol(&(solv[1]));
   }
   /* for structure */
   if (genprob.probtyp == prb_structure)
   {
      if (genprob.numfld!=1) dserror("numfld != 1 for Structural Problem");

      solv = (SOLVAR*)CCACALLOC(genprob.numfld,sizeof(SOLVAR));

      solv[genprob.numsf].fieldtyp = structure;
      inpctrsol(&(solv[genprob.numsf]));
   }
   /* for optimisation */
   if (genprob.probtyp == prb_opt)
   {
      if (genprob.numfld!=1) dserror("numfld != 1 for Structural Problem");

      solv = (SOLVAR*)CCACALLOC(genprob.numfld,sizeof(SOLVAR));

      solv[0].fieldtyp = structure;
      inpctrsol(&(solv[0]));
   }
   /* for fluid */
   if (genprob.probtyp == prb_fluid)
   {
      solv = (SOLVAR*)CCACALLOC(genprob.numfld,sizeof(SOLVAR));

      solv[genprob.numff].fieldtyp = fluid;
      inpctrsol(&(solv[genprob.numff]));

      if (genprob.numfld==2)
      {
        solv[genprob.numaf].fieldtyp = ale;
	inpctrsol(&(solv[genprob.numaf]));
      }
   }
   /* for (plain) ALE */
   if (genprob.probtyp == prb_ale)
   {
      if (genprob.numfld!=1) dserror("numfld != 1 for Ale Problem");

      solv = (SOLVAR*)CCACALLOC(genprob.numfld,sizeof(SOLVAR));

      solv[genprob.numaf].fieldtyp = ale;
      inpctrsol(&(solv[genprob.numaf]));
   }
#ifdef D_TSI
   /* for TSI */
   if (genprob.probtyp == prb_tsi)
   {
      if (genprob.numfld != 2)
      {
        dserror("numfld != 2 for TSI: impossible");
      }

      /* allocate memory for solver entries */
      solv = (SOLVAR *) CCACALLOC(genprob.numfld, sizeof(SOLVAR));

      /* set solver details of structure solver */
      solv[genprob.numsf].fieldtyp = structure;
      inpctrsol(&(solv[genprob.numsf]));

      /* set solver details of thermal solver */
      solv[genprob.numtf].fieldtyp = thermal;
      inpctrsol(&(solv[genprob.numtf]));
   }
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inpctr */




/*----------------------------------------------------------------------*
 | input of general problem data                          m.gee 2/01    |
 *----------------------------------------------------------------------*/
void inpctrprob()
{
INT  ierr;
INT  flag = 0;
INT  restart;
char buffer[50];
#ifdef DEBUG
dstrc_enter("inpctrprob");
#endif

/*------------------------------------------ set initial field numbers */
genprob.numsf=-1;  /* structure fields */
genprob.numff=-1;  /* fluid fields */
genprob.numaf=-1;  /* ALE fields */
genprob.numls=-1;
#ifdef D_TSI
genprob.numtf=-1;  /* thermal fields */
#endif

if (frfind("-PROBLEM SIZE")==0)
{
  dserror("frfind: PROBLEM SIZE not in input file");
}
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
  frint("ELEMENTS", &(genprob.nele),&ierr);
  frint("NODES",    &(genprob.nnode),&ierr);
  frint("DIM",      &(genprob.ndim),&ierr);
  frint("MATERIALS",&(genprob.nmat),&ierr);
  frint("NUMDF",    &(genprob.numdf),&ierr);

  frread();
}
/*------------------------------------------------------ check values */
if (genprob.nmat<=0)
   dserror("No Material defined!");
/*-------------------------------------- default value for monitoring */
ioflags.monitor=0;
/*----------------------------------------- default values for output */
ioflags.output_out          =0;
ioflags.output_gid          =0;
ioflags.output_bin          =0;
ioflags.struct_disp         =0;
ioflags.struct_stress       =0;
ioflags.struct_stress_smo   =0;
ioflags.struct_sm_disp      =0;
ioflags.struct_sm_stress    =0;
ioflags.fluid_sol           =0;
ioflags.fluid_stress        =0;
ioflags.fluid_vis           =0;
ioflags.ale_disp            =0;
#ifdef D_TSI
ioflags.therm_temper        =0;
ioflags.therm_heatflux      =0;
#endif
ioflags.relative_displ      =0;
ioflags.output_dis          =0;


/* set default values */
genprob.usetrilinosalgebra=0;
/* read values */
if (frfind("-PROBLEM TYP")==0) dserror("frfind: PROBLEM TYP not in input file");
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
  /* read problem type */
  frchar("PROBLEMTYP",buffer,   &ierr);
  if (ierr==1)
  {
    if (frwordcmp("Structure"                  ,buffer)==0) genprob.probtyp = prb_structure;
    if (frwordcmp("Fluid"                      ,buffer)==0) genprob.probtyp = prb_fluid;
    if (frwordcmp("Fluid_Structure_Interaction",buffer)==0) genprob.probtyp = prb_fsi;
    if (frwordcmp("Structure_Structure_Interaction",buffer)==0) genprob.probtyp = prb_ssi;
    if (frwordcmp("Optimisation"               ,buffer)==0) genprob.probtyp = prb_opt;
    if (frwordcmp("Ale"                        ,buffer)==0) genprob.probtyp = prb_ale;
#ifdef D_TSI
    if (frwordcmp("Thermal_Structure_Interaction",buffer)==0) genprob.probtyp = prb_tsi;
#endif
  }

  /* read time type : time-dependent/time-independent analysis */
  frchar("TIMETYP"   ,buffer,            &ierr);
  if (ierr==1)
  {
    if (frwordcmp("Static" ,buffer)==0) genprob.timetyp=time_static;
    if (frwordcmp("Dynamic",buffer)==0) genprob.timetyp=time_dynamic;
  }

  /* restarting flag */
  frint("RESTART"    ,&(restart),&ierr);
  if (ierr==1)
  {
     if (genprob.restart==0 && restart != 0)
     {
       dserror("Restart defined in input file but not as program argument!\n");
     }
     genprob.restart=restart;
  }

  /* read number of fields */
  frint("NUMFIELD",&(genprob.numfld),&ierr);

  frint("GRADERW",&(genprob.graderw),&ierr);

  frint("MULTISC_STRUCT",&(genprob.multisc_struct),&ierr);

  frchar("ALGEBRA",buffer,   &ierr);
  if (ierr)
  {
    if (frwordcmp("Trilinos",buffer)==0) genprob.usetrilinosalgebra = 1;
    if (frwordcmp("ccarat"  ,buffer)==0) genprob.usetrilinosalgebra = 0; /* default */
#ifndef TRILINOS_PACKAGE
    if (genprob.usetrilinosalgebra) dserror("ALGEBRA Trilinos in input file requires -DTRILINOS_PACKAGE");
#endif

  }

  frread();
}

/*----------------- set field numbers depending on problem type and numfld */
if (genprob.probtyp==prb_fsi)
{
   genprob.numsf=0;
   genprob.numff=1;
   genprob.numaf=2;
}
if (genprob.probtyp==prb_fluid)
{
   genprob.numff=0;
   if (genprob.numfld==2) genprob.numaf=1;
}
if (genprob.probtyp==prb_ale) genprob.numaf=0;
if (genprob.probtyp==prb_structure) genprob.numsf=0;
#ifdef D_SSI
if (genprob.probtyp==prb_ssi)
{
   genprob.numsf=0;
}
#endif
#ifdef D_TSI
if (genprob.probtyp == prb_tsi)
{
   genprob.numsf = 0;  /* structural field index */
   genprob.numtf = 1;  /* thermal field index */
}
#endif

/*----------------------------------------------------------------------*/
/* input / output file choices */
if (frfind("---IO")==1)
{
  frread();
  ioflags.steps_per_file = 1000;
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    frreadyes("OUTPUT_OUT",&ioflags.output_out);
    frreadyes("OUTPUT_GID",&ioflags.output_gid);
    frreadyes("OUTPUT_BIN",&ioflags.output_bin);
    frreadyes("STRUCT_DISP",&ioflags.struct_disp);
    frreadyes("STRUCT_STRESS",&ioflags.struct_stress);
    frreadyes("STRUCT_STRESS_SMO",&ioflags.struct_stress_smo);
    frreadyes("STRUCT_SM_DISP",&ioflags.struct_sm_disp);
    frreadyes("STRUCT_SM_STRESS",&ioflags.struct_sm_stress);
    frreadyes("FLUID_SOL",&ioflags.fluid_sol);
    frreadyes("FLUID_STRESS",&ioflags.fluid_stress);
    frreadyes("FLUID_VIS",&ioflags.fluid_vis);
    frreadyes("ALE_DISP",&ioflags.ale_disp);
#ifdef D_TSI
    frreadyes("THERM_TEMPERATURE",&ioflags.therm_temper);
    frreadyes("THERM_HEATFLUX",&ioflags.therm_heatflux);
#endif

    /* the old keywords (only for backwards compatibility) */
    ierr = frreadyes("STRUCT_DISP_FILE",&flag);
    if (ierr && flag)
    {
      ioflags.output_out=1;
      ioflags.struct_disp=1;
    }

    ierr = frreadyes("STRUCT_STRESS_FILE",&flag);
    if (ierr && flag)
    {
      ioflags.output_out=1;
      ioflags.struct_stress=1;
    }

    ierr = frreadyes("STRUCT_DISP_GID",&flag);
    if (ierr && flag)
    {
      ioflags.output_gid=1;
      ioflags.struct_disp=1;
    }

    ierr = frreadyes("STRUCT_STRESS_GID",&flag);
    if (ierr && flag)
    {
      ioflags.output_gid=1;
      ioflags.struct_stress=1;
    }

    ierr = frreadyes("STRUCT_STRESS_GID_SMO",&flag);
    if (ierr && flag)
    {
      ioflags.output_gid=1;
      ioflags.struct_stress_smo=1;
    }

    ierr = frreadyes("STRUCT_SM_DISP_GID",&flag);
    if (ierr && flag)
    {
      ioflags.output_gid=1;
      ioflags.struct_sm_disp=1;
    }

    ierr = frreadyes("STRUCT_SM_STRESS_GID",&flag);
    if (ierr && flag)
    {
      ioflags.output_gid=1;
      ioflags.struct_sm_stress=1;
    }

    ierr = frreadyes("FLUID_STRESS_GID",&flag);
    if (ierr && flag)
    {
      ioflags.output_gid=1;
      ioflags.fluid_stress=1;
    }

    ierr = frreadyes("FLUID_SOL_GID",&flag);
    if (ierr && flag)
    {
      ioflags.output_gid=1;
      ioflags.fluid_sol=1;
    }

    ierr = frreadyes("FLUID_SOL_FILE",&flag);
    if (ierr && flag)
    {
      ioflags.output_out=1;
      ioflags.fluid_sol=1;
    }

    frreadyes("FLUID_VIS_FILE",&ioflags.fluid_vis);

    ierr = frreadyes("ALE_DISP_FILE",&flag);
    if (ierr && flag)
    {
      ioflags.output_out=1;
      ioflags.ale_disp=1;
    }

    ierr = frreadyes("ALE_DISP_GID",&flag);
    if (ierr && flag)
    {
      ioflags.output_gid=1;
      ioflags.ale_disp=1;
    }
    /* End of old keywords */



    frint("RELATIVE_DISP_NUM",&(ioflags.relative_displ),&ierr);
    frint("FILESTEPS",&(ioflags.steps_per_file),&ierr);
    frread();
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inpctrprob */




/*----------------------------------------------------------------------*
 | input of static problem data                           m.gee 2/01    |
 *----------------------------------------------------------------------*/
void inpctrstat()
{
INT  ierr;
char buffer[50];
INT  counter;
#ifdef DEBUG
dstrc_enter("inpctrstat");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------- allocate a structure STATIC */
statvar = (STATIC_VAR*)CCACALLOC(1,sizeof(STATIC_VAR));
if (!statvar) dserror("Allocation of STATIC failed");
/*----------------------------- default: results written every loadstep */
statvar->resevry_disp=1;
statvar->resevry_stress=1;
statvar->resevery_restart=1;
statvar->isrelstepsize=0;
statvar->multiscale=0;
/*------------------------------------------------------- start reading */
if (frfind("-STATIC")==1)
{
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    frchk("LINEAR",&ierr);
    if (ierr==1) {statvar->linear=1;statvar->nonlinear=0;}
    frchk("NONLINEAR",&ierr);
    if (ierr==1) {statvar->nonlinear=1;statvar->linear=0;}
    /*------------------------- read if multiscale model */
    frchar("MULTISCALE",buffer,&ierr);
    if (ierr)
    {
      if (frcheckyes(buffer)) statvar->multiscale=1;
    }
    /*------------------------- read for typ of kinematic sh 03/03 */
    frchar("KINTYP",buffer,&ierr);
    if (ierr==1)
    {
      if (frwordcmp(buffer,"Geo_Lin")==0)
        statvar->kintyp = geo_linear;
      if (frwordcmp(buffer,"Upd_Lagr")==0)
        statvar->kintyp = upd_lagr;
      if (frwordcmp(buffer,"Tot_Lagr")==0)
        statvar->kintyp = tot_lagr;
    }
    /*------------------------- read for typ of pathfollowing technique */
    frchar("NEWTONRAPHSO",buffer,&ierr);
    if (ierr==1)
    {
      if (frwordcmp(buffer,"Displacement_Control")==0)
        statvar->nr_controltyp = control_disp;
      if (frwordcmp(buffer,"Load_Control")==0)
        statvar->nr_controltyp = control_load;
      if (frwordcmp(buffer,"Arc_Control")==0)
        statvar->nr_controltyp = control_arc;
      if (frwordcmp(buffer,"none")==0)
        statvar->nr_controltyp = control_none;
    }
    /*---------------------------------------- read for arcscaling flag */
    frreadyes("IARC", &(statvar->iarc));
    /*------------------------------ read for sign-changing-by-csp flag */
    frreadyes("SIGNCHCSP", &(statvar->signchcsp));
    /*-------------------------------------------------- read variables */
    frint("NUMSTEP",&(statvar->nstep)  ,&ierr);
    frint("MAXITER",&(statvar->maxiter),&ierr);
    frint("RESEVRYDISP",&(statvar->resevry_disp)  ,&ierr);
    frint("RESEVRYSTRS",&(statvar->resevry_stress),&ierr);
    frint("RESTARTEVRY",&(statvar->resevery_restart),&ierr);
    frint("GRADERW",&(statvar->graderw),&ierr);

    frdouble("TOLRESID",&(statvar->tolresid),&ierr);
    frdouble("TOLDISP" ,&(statvar->toldisp) ,&ierr);
    frdouble("STEPSIZE",&(statvar->stepsize),&ierr);
    frdouble("ARCSCL"  ,&(statvar->arcscl)  ,&ierr);
    /*-------------------------------------------------------------------*/
    frreadyes("VARSTEPSI", &(statvar->isrelstepsize));
    /*-------------------------------------------------------------------*/
    frread();
  }
}
/*--------------------------------- in nonlinear case find control node */
if (statvar->nonlinear)
{
  /*---------------------------------------------------- start reading */
  if (frfind("-CONTROL NODE")==1)
  {
    frread();
    while(strncmp(allfiles.actplace,"------",6)!=0)
    {
      frint("NODE",&(statvar->control_node_global),&ierr);
      if (ierr) statvar->control_node_global--;
      frint("DOF",&(statvar->control_dof),&ierr);
      frread();
    }
  }
  frrewind();
}
/*----------------------------------------------------------------*/
if (ioflags.relative_displ>0)
{
  if (frfind("-RELATIVE DISPLACEMENT NODES")==1)
  {
    frread();
    counter=0;
    while(strncmp(allfiles.actplace,"------",6)!=0)
    {
      frint("NODE",&(statvar->reldisnode_ID[counter]),&ierr);
      if (ierr) statvar->reldisnode_ID[counter]--;
      frint("DOF",&(statvar->reldis_dof[counter]),&ierr);
      counter++;
      frread();
    }
  }
}
/*----------------------------------------------------------------*/
if (statvar->isrelstepsize==1)
{
  if (frfind("-VARIABLE STEP SIZES")==1)
  {
    frread();
    counter=0;
    while(strncmp(allfiles.actplace,"------",6)!=0)
    {
      frint("NSTEP",&(statvar->actstep[counter]),&ierr);
      frdouble("ACTSTEPSIZE",&(statvar->actstepsize[counter]),&ierr);
      counter++;
      frread();
    }
    statvar->numcurve = counter;
    if(counter>19)  dserror("not more than 20 different stepsizes!");
  }
}
/*----------------------------------------------------------------*/
#ifdef D_MLSTRUCT
if (genprob.multisc_struct == 1)
{
  if (frfind("--SMVALUES")==1)
  {
     frread();
     while(strncmp(allfiles.actplace,"------",6)!=0)
     {
        frdouble("EPS_EQUIVAL",&(statvar->eps_equiv),&ierr);
        frread();
     }
  }
}
#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inpctrstat */


/*----------------------------------------------------------------------*
 | input of data for eigensolution                           al 8/02    |
 *----------------------------------------------------------------------*/
void inpctreig()
{
INT    i;
FIELD *actfield;
#ifdef DEBUG
dstrc_enter("inpctreig");
#endif
/*----------------------------------------------------------------------*/
alleig = (ALLEIG*)CCACALLOC(1,sizeof(ALLEIG));
if (!alleig) dserror("Allocation of ALLEIG failed");
/*----------------------------------------------------------------------*/
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
   switch(actfield->fieldtyp)
   {
   case fluid:
   break;
   case ale:
   break;
   case structure:
     inpctr_eig_struct(alleig);
   break;
#ifdef D_TSI
   case thermal:
   break;
#endif
   case none:
      dserror("Cannot find type of field");
   break;
   default:
      dserror("Cannot find type of field");
   break;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inpctreig */

/*----------------------------------------------------------------------*
 | input of eigensolution problem data                       al 8/02    |
 *----------------------------------------------------------------------*/
void inpctr_eig_struct(ALLEIG *alleig)
{

INT    ierr;
char   buffer[50];
#ifdef DEBUG
dstrc_enter("inpctr_eig_struct");
#endif

if (frfind("--EIGENVALUE ANALYSIS")==0) goto end;
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("SOLTYP"   ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"SUBSPACE")==0)
          alleig->soltyp=subspace;
      else
          alleig->soltyp=eig_none;
   }
/*--------------read INT */
   frint("STURM ",&(alleig->sturm ),&ierr);
   frint("SUBTYP",&(alleig->subtyp),&ierr);
   frint("IFSH"  ,&(alleig->ifsh  ),&ierr);
   frint("ILMP"  ,&(alleig->ilmp  ),&ierr);
   frint("RANGE" ,&(alleig->range ),&ierr);
   frint("NUMVEC",&(alleig->numvec),&ierr);
   frint("NROOT" ,&(alleig->nroot ),&ierr);
   frint("ITEMAX",&(alleig->itemax),&ierr);
   frint("IFCTR" ,&(alleig->ifctr ),&ierr);
/*--------------read DOUBLE */
   frdouble("TOLEIG",&(alleig->toleig),&ierr);
   frdouble("SHIFT" ,&(alleig->shift ),&ierr);
   frdouble("BOULO" ,&(alleig->boulo ),&ierr);
   frdouble("BOUUP" ,&(alleig->bouup ),&ierr);
   frread();
}
frrewind();

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inpctr_eig_struct */


/*----------------------------------------------------------------------*
 | input of dynamic problem data                          m.gee 2/01    |
 |                                                         genk 3/01    |
 *----------------------------------------------------------------------*/
void inpctrdyn()
{
INT    i;
FIELD *actfield;
#ifdef DEBUG
dstrc_enter("inpctrdyn");
#endif
/*----------------------------------------------------------------------*/
if (genprob.probtyp==prb_fsi ||
   (genprob.probtyp==prb_fluid &&  genprob.numfld>1))
   alldyn = (ALLDYNA*)CCACALLOC((genprob.numfld+1),sizeof(ALLDYNA));
else if (genprob.probtyp==prb_ssi)
   alldyn = (ALLDYNA*)CCACALLOC((genprob.numfld+1),sizeof(ALLDYNA));
#ifdef D_TSI
else if (genprob.probtyp == prb_tsi)
   alldyn = (ALLDYNA *) CCACALLOC((genprob.numfld+1), sizeof(ALLDYNA));
#endif
else
   alldyn = (ALLDYNA*)CCACALLOC(genprob.numfld,sizeof(ALLDYNA));
if (!alldyn) dserror("Allocation of ALLDYNA failed");
/*----------------------------------------------------------------------*/
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
   switch(actfield->fieldtyp)
   {
       case fluid:
#ifdef D_FLUID
         alldyn[i].fdyn = (FLUID_DYNAMIC*)CCACALLOC(1,sizeof(FLUID_DYNAMIC));
         inpctr_dyn_fluid(alldyn[i].fdyn);
#else
         dserror("General FLUID problem not defined in Makefile!!!");
#endif
         break;
       case ale:
#ifdef D_ALE
         alldyn[i].adyn = (ALE_DYNAMIC*)CCACALLOC(1,sizeof(ALE_DYNAMIC));
         inpctr_dyn_ale(alldyn[i].adyn);
#else
         dserror("General ALE problem not defined in Makefile!!!");
#endif
         break;
       case structure:
         alldyn[i].sdyn = (STRUCT_DYNAMIC*)CCACALLOC(1,sizeof(STRUCT_DYNAMIC));
         inpctr_dyn_struct(alldyn[i].sdyn);
         break;
#ifdef D_TSI
       case thermal:
         alldyn[i].tdyn = (THERM_DYNAMIC*)CCACALLOC(1,sizeof(THERM_DYNAMIC));
         /* this is empty right now, come back later */
         break;
#endif
       case none:
         dserror("Cannot find type of field");
         break;
       default:
         dserror("Cannot find type of field");
         break;
   }
}
/*----------------------------------------------------------------------*/
 if (genprob.probtyp==prb_fsi ||
     (genprob.probtyp==prb_fluid &&  genprob.numfld>1))
 {
#ifdef D_FSI
   alldyn[i].fsidyn = (FSI_DYNAMIC*)CCACALLOC(1,sizeof(FSI_DYNAMIC));
   inpctr_dyn_fsi(alldyn[i].fsidyn);
#else
   dserror("General FSI problem not defined in Makefile!!!");
#endif
 }
/*----------------------------------------------------------------------*/
if (genprob.probtyp==prb_ssi)
{
#ifdef D_SSI
   alldyn[i].ssidyn = (SSI_DYNAMIC*)CCACALLOC(1,sizeof(SSI_DYNAMIC));
   inpctr_dyn_ssi(alldyn[i].ssidyn);
#else
   dserror("General SSI problem not defined in Makefile!!!");
#endif
}

/*----------------------------------------------------------------------*/
if (genprob.probtyp == prb_tsi)
{
#ifdef D_TSI
  /* We allocated an additional ALLDYNA component, the (numfld+1)th,
     this is now set with the interaction dyn controls */
  /* counter i==genprob.numfld */
  alldyn[i].tsidyn = (TSI_DYNAMIC*) CCACALLOC(1, sizeof(TSI_DYNAMIC));
  inpctr_dyn_tsi(alldyn[i].tsidyn);
#else
  dserror("General TSI problem not defined in Makefile!!!");
#endif
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
 dstrc_exit();
#endif
 return;
} /* end of inpctrdyn */



/*----------------------------------------------------------------------*
 | input of dynamic problem data  for field structure     m.gee 2/02    |
 *----------------------------------------------------------------------*/
void inpctr_dyn_struct(STRUCT_DYNAMIC *sdyn)
{
INT    ierr;
char   buffer[50];
#ifdef DEBUG
dstrc_enter("inpctr_dyn_struct");
#endif

sdyn->updevry_disp=1;
sdyn->updevry_stress=1;
sdyn->res_write_evry=1;
sdyn->eigen=0;

if (frfind("-STRUCTURAL DYNAMIC")==0) goto end;
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("DYNAMICTYP",buffer,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"Centr_Diff")==0)
      {
        sdyn->Typ = centr_diff;
      }
      else if (frwordcmp(buffer,"Gen_Alfa")==0)
      {
        sdyn->Typ = gen_alfa;
      }
      else if (frwordcmp(buffer,"Gen_EMM")==0)
      {
        sdyn->Typ = Gen_EMM;
      }
      else
      {
        sdyn->Typ = gen_alfa;
        printf("DYNAMICTYP unknown (the time integration scheme) "
          "using default: Generalised-alpha");
      }
   }
   frreadyes("DAMPING",&(sdyn->damp));
   frchar("ITERATION"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"full")==0 ||
          frwordcmp(buffer,"Full")==0 ||
          frwordcmp(buffer,"FULL")==0 )

          sdyn->iter=1;
      else
          sdyn->iter=0;
   }

/*-----------------read INT */
   frint("EIGEN"      ,&(sdyn->eigen)         ,&ierr);
   frint("NUMSTEP"    ,&(sdyn->nstep)         ,&ierr);
   frint("MAXITER"    ,&(sdyn->maxiter)       ,&ierr);
   frint("RESEVRYDISP",&(sdyn->updevry_disp)  ,&ierr);
   frint("RESEVRYSTRS",&(sdyn->updevry_stress),&ierr);
   frint("RESTARTEVRY",&(sdyn->res_write_evry),&ierr);
   frint("CONTACT"    ,&(sdyn->contact)       ,&ierr);
#ifdef WALLCONTACT
   frint("CET_flag",   &(contact.CET_flag)    ,&ierr);
   frint("FR_flag",    &(contact.FR_flag)     ,&ierr);
#endif
/*--------------read DOUBLE */
   frdouble("TIMESTEP",&(sdyn->dt)     ,&ierr);
#ifdef WALLCONTACT
   if (ierr) contact.dt = sdyn->dt;
#endif
   frdouble("MAXTIME" ,&(sdyn->maxtime),&ierr);
   frdouble("BETA"    ,&(sdyn->beta)   ,&ierr);
   frdouble("GAMMA"   ,&(sdyn->gamma)  ,&ierr);
   frdouble("ALPHA_M" ,&(sdyn->alpha_m),&ierr);
   frdouble("ALPHA_F" ,&(sdyn->alpha_f),&ierr);
#ifdef GEMM
   frdouble("XSI"     ,&(sdyn->xsi)    ,&ierr);
#endif
#ifdef WALLCONTACT
   frdouble("NPP"     ,&(contact.n_pen_par) , &ierr);
   frdouble("TPP"     ,&(contact.t_pen_par) , &ierr);
   frdouble("FR_COEF" ,&(contact.fr_coef)   , &ierr);
#endif
   frdouble("M_DAMP"  ,&(sdyn->m_damp) ,&ierr);
   frdouble("K_DAMP"  ,&(sdyn->k_damp) ,&ierr);
   frdouble("TOLDISP" ,&(sdyn->toldisp),&ierr);

/*------ read time adaption */
   frint("TIMEADAPT"  ,&(sdyn->timeadapt),&ierr);
   frint("ITWANT"     ,&(sdyn->itwant)   ,&ierr);
   frdouble("MAXDT"   ,&(sdyn->maxdt)    ,&ierr);
   frdouble("RESULTDT",&(sdyn->resultdt) ,&ierr);

   frread();
}
frrewind();
/*----------------------------------------------------------------------*/

end:
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inpctr_dyn_struct */

#ifdef D_FLUID

/*!---------------------------------------------------------------------
\brief input of the FLUID DYNAMIC block in the input-file

<pre>                                                         genk 03/02

In this routine the data in the FLUID DYNAMIC block of the input file
are read and stored in fdyn

</pre>
\param  *fdyn 	  FLUID_DATA       (o)
\return void

------------------------------------------------------------------------*/
void inpctr_dyn_fluid(FLUID_DYNAMIC *fdyn)
{
INT    ierr;
INT    i;
INT    thetafound=0;
char   buffer[50];
#ifdef DEBUG
dstrc_enter("inpctr_dyn_fluid");
#endif

/*-------------------------------------------------- set default values */
/* general control variables of fluid dynamics */
fdyn->dyntyp=dyntyp_nln_time_int; /* default: nonlin. time integration */
fdyn->init=0;
fdyn->iop=2;            /* default: one-step-theta */
fdyn->iops=2;           /* default: one-step-theta */
fdyn->mlfem = 0;	/* default: no multi-level algorithm */
fdyn->numdf=genprob.ndim+1;
fdyn->ite=1;            /* default: fixed-point like */
fdyn->itnorm=fncc_L2;   /* default: L2-norm */
fdyn->itemax=10;
fdyn->stchk=-1;         /* default: no steady state check */
fdyn->stnorm=fnst_L2;   /* default: L2-norm */
/* output flags */
fdyn->uppss=1;
fdyn->upout=1;
fdyn->upres=1;
fdyn->uprestart=20;
fdyn->resstep=0;
/* time stepping flags and variables */
fdyn->itnum=0;
fdyn->nstep=1;
fdyn->step=0;
fdyn->nums=0;
/* time integration variables */
fdyn->maxtime=1000.0;
fdyn->dt=0.01;
fdyn->max_dt = 1.0;	/* default: maximal time step size to 1.0 */
fdyn->min_dt = 0.0;	/* default: minimal time step size to 0.0 */
fdyn->theta=0.66;
fdyn->thetas=ONE;
fdyn->alpha_f=ONE;
fdyn->alpha_m=ONE;
fdyn->lte = EPS3;	/* default: local truncation error for adapt. time*/
/* special facilities flags */
fdyn->iprerhs=1;        /* this value is not in the input file! */
fdyn->viscstr=0;        /* default: do not include viscose stresses */
fdyn->freesurf=0;       /* default: no free surface */
fdyn->surftens=0;       /* default: do not include surface tension */
fdyn->checkarea=0;
fdyn->liftdrag=ld_none;
fdyn->adaptive=0;	/* default: no adaptive time stepping */
fdyn->stresspro=0;      /* default: do no stress projection step */

/* turbulence flags */
fdyn->turbu=0;
fdyn->dis_capt=0;
fdyn->itemax_ke=100;
fdyn->stepke=0;


/* tolerances */
fdyn->sttol=EPS6;
fdyn->ittol=EPS6;


if (frfind("-FLUID DYNAMIC")==0) goto end;
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("DYNAMICTYP",buffer,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"Projection_Method")==0)
	fdyn->dyntyp=dyntyp_pm_discont;
      else if (frwordcmp(buffer,"PM_Laplace")==0)
	fdyn->dyntyp=dyntyp_pm_cont_laplace;
      else if (frwordcmp(buffer,"Nlin_Time_Int")==0)
	fdyn->dyntyp=dyntyp_nln_time_int;
      else
	dserror("DYNAMICTYP unknown");
   }
   frchar("TIMEINTEGR"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"Stationary")==0)
         fdyn->iop=0;
      else if (frwordcmp(buffer,"Gen_Alfa")==0)
         fdyn->iop=1;
      else if (frwordcmp(buffer,"Gen_Alpha")==0)
         fdyn->iop=1;
      else if (frwordcmp(buffer,"One_Step_Theta")==0)
         fdyn->iop=4;
      else if (frwordcmp(buffer,"BDF2")==0)
         fdyn->iop=7;
      else
         dserror("TIMEINTEGR-Type unknown");
   }
   frchar("STARTINGALGO"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"One_Step_Theta")==0)
         fdyn->iops=4;
      else
         dserror("STARTINGALGO unknown");
   }
   frchar("NONLINITER"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"fixed_point_like")==0)
         fdyn->ite=1;
      else if (frwordcmp(buffer,"Newton")==0)
         fdyn->ite=2;
      else
         dserror("NONLINITER-Type unknown!");
   }
   frchar("CONVCHECK"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frcheckno(buffer))
	 fdyn->itnorm=fncc_no;
      else if (frwordcmp(buffer,"L_infinity_norm")==0)
         fdyn->itnorm=fncc_Linf;
      else if (frwordcmp(buffer,"L_1_norm")==0)
         fdyn->itnorm=fncc_L1;
      else if (frwordcmp(buffer,"L_2_norm")==0)
         fdyn->itnorm=fncc_L2;
      else
         dserror("Norm for CONVCHECK unknown");
   }
   frchar("STEADYCHECK"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frcheckno(buffer))
	 fdyn->stnorm=fnst_no;
      else if (frwordcmp(buffer,"L_infinity_norm")==0)
         fdyn->stnorm=fnst_Linf;
      else if (frwordcmp(buffer,"L_1_norm")==0)
         fdyn->stnorm=fnst_L1;
      else if (frwordcmp(buffer,"L_2_norm")==0)
         fdyn->stnorm=fnst_L2;
      else
         dserror("Norm for STEADYCHECK unknown");
   }
   frchar("INITIALFIELD"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"zero_field")==0)
         fdyn->init=0;
      else if (frwordcmp(buffer,"field_from_file")==0)
         fdyn->init=1;
      else if (frwordcmp(buffer,"field_by_function")==0)
         fdyn->init=-1;
      else if (frwordcmp(buffer,"SOLWAVE")==0)
         fdyn->init=6;
      else if (frwordcmp(buffer,"WAVEBREAKING")==0)
         fdyn->init=7;
      else if (frwordcmp(buffer,"BELTRAMI-FLOW")==0)
         fdyn->init=8;
      else if (frwordcmp(buffer,"KIM-MOIN-FLOW")==0)
         fdyn->init=9;
      else if (frwordcmp(buffer,"BREAKING-DAM")==0)
         fdyn->init=10;
      else
         dserror("INITIALFIELD unknown!");
   }
   frreadyes("VISCSTRESS",&(fdyn->viscstr));
   frreadyes("STRESSPRO",&(fdyn->stresspro));
   frchar("FREESURFACE"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frcheckno(buffer))
         fdyn->freesurf=0;
      else if (frwordcmp(buffer,"loclag_exp")==0)
         fdyn->freesurf=1;
      else if (frwordcmp(buffer,"loclag_imp")==0)
         fdyn->freesurf=2;
      else if (frwordcmp(buffer,"hf_vert_sep")==0)
         fdyn->freesurf=3;
      else if (frwordcmp(buffer,"hf_vert_imp")==0)
         fdyn->freesurf=5;
      else if (frwordcmp(buffer,"genfs")==0)
         fdyn->freesurf=6;
      else
         dserror("Parameter FREESURFACE unknown!\n");
   }
   frreadyes("SURFTENSION",&(fdyn->surftens));
   frreadyes("CHECKAREA",&(fdyn->checkarea));
   frchar("LIFTDRAG"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frcheckyes(buffer))
         fdyn->liftdrag=ld_stress;
      else if (frwordcmp(buffer,"nodeforce")==0 ||
               frwordcmp(buffer,"NODEFORCE")==0 ||
	       frwordcmp(buffer,"Nodeforce")==0    )
         fdyn->liftdrag=ld_nodeforce;
      else if (frwordcmp(buffer,"stress")==0 ||
               frwordcmp(buffer,"STRESS")==0 ||
	       frwordcmp(buffer,"Stress")==0    )
         fdyn->liftdrag=ld_stress;
      else if (frcheckno(buffer))
         fdyn->liftdrag=ld_none;
      else
         dserror("LIFTDRAG unknown!\n");
   }
   frchar("TURBULENCE"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frcheckno(buffer))
         fdyn->turbu=0;
      else if (frwordcmp(buffer,"algebraic")==0)
         fdyn->turbu=1;
      else if (frwordcmp(buffer,"kappa-eps")==0)
         fdyn->turbu=2;
      else if (frwordcmp(buffer,"kappa-omega")==0)
         fdyn->turbu=3;
      else
         dserror("TURBULENCE unknown!");
   }
   frreadyes("DISC_CAPT",&(fdyn->dis_capt));
   frreadyes("ADAPT_TIME",&(fdyn->adaptive));
   frchar("CONVECTERM"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"convective")==0)
         fdyn->conte=0;
      else if (frwordcmp(buffer,"divergence")==0)
         fdyn->conte=1;
      else if (frwordcmp(buffer,"skew_symmetric")==0)
         fdyn->conte=2;
      else
         dserror("Unknown convective term!");
   }
   frchar("VISCTERM"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"conventional")==0)
         fdyn->vite=0;
      else if (frwordcmp(buffer,"stress_divergence")==0)
         fdyn->vite=1;
      else
         dserror("Unknown viscous term!");
   }
   frchar("SUBGRIDVISC"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frcheckno(buffer))
         fdyn->sgvisc=0;
      else if (frwordcmp(buffer,"artificial")==0)
         fdyn->sgvisc=1;
      else if (frwordcmp(buffer,"Smagorinsky")==0)
         fdyn->sgvisc=2;
      else
         dserror("Unknown subgrid viscosity!");
   }

/*--------------read INT */
   frint("UPPSS"      ,&(fdyn->uppss)     ,&ierr);
   frint("UPOUT"      ,&(fdyn->upout)     ,&ierr);
   frint("UPRES"      ,&(fdyn->upres)     ,&ierr);
   frint("RESSTEP"    ,&(fdyn->resstep)   ,&ierr);
   frint("RESTARTEVRY",&(fdyn->uprestart) ,&ierr);
   frint("NUMSTEP"    ,&(fdyn->nstep)     ,&ierr);
   frint("STEADYSTEP" ,&(fdyn->stchk)     ,&ierr);
   frint("NUMSTASTEPS" ,&(fdyn->nums)     ,&ierr);
   frint("STARTFUNCNO" ,&i                ,&ierr);
   if (ierr==1)
   {
      if (fdyn->init==-1)
          fdyn->init=i;
   }
   frint("ITEMAX" ,&(fdyn->itemax) ,&ierr);

/*--------------read DOUBLE */
   frdouble("TIMESTEP" ,&(fdyn->dt)     ,&ierr);
   frdouble("MAXTIME"  ,&(fdyn->maxtime),&ierr);
   frdouble("ALPHA_M"  ,&(fdyn->alpha_m),&ierr);
   frdouble("ALPHA_F"  ,&(fdyn->alpha_f),&ierr);
   if (thetafound==0)
   {
      frdouble("THETA"    ,&(fdyn->theta)  ,&ierr);
      if (ierr==1) thetafound++;
   }
   frdouble("CONVTOL"     ,&(fdyn->ittol)  ,&ierr);
   frdouble("STEADYTOL"   ,&(fdyn->sttol) ,&ierr);
   frdouble("START_THETA" ,&(fdyn->thetas),&ierr);
   frdouble("INT_LENGHT"  ,&(fdyn->lenght),&ierr);
   frdouble("ROUGHTNESS"  ,&(fdyn->rought),&ierr);
   frdouble("SC_COORD_X"  ,&(fdyn->coord_scale[0]),&ierr);
   frdouble("SC_COORD_Y"  ,&(fdyn->coord_scale[1]),&ierr);
   frdouble("MAX_DT"      ,&(fdyn->max_dt),&ierr);
   frdouble("MIN_DT"      ,&(fdyn->min_dt),&ierr);
   frdouble("LOC_TRUN_ERR",&(fdyn->lte)   ,&ierr);
   frdouble("SMAGCONST",&(fdyn->smagcon),&ierr);

   frread();
}
frrewind();

/*----------------------------------------------------------------------*/
if (frfind("-MULTILEVEL FLUID DYNAMIC")==0) goto end;
frread();
/* allocate fluid mulitlevel variables */
fdyn->mlvar = (FLUID_DYN_ML*)CCACALLOC(1,sizeof(FLUID_DYN_ML));

while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("MULEVELFEM"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no",2)==0)
         fdyn->mlfem=0;
      else if (strncmp(buffer,"yes",3)==0)
         fdyn->mlfem=1;
      else
         dserror("unknown Multilevel");
   }
   frchar("TRANSTERM"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no_approx",9)==0)
         fdyn->mlvar->transterm=0;
      else if (strncmp(buffer,"neglect_sst",11)==0)
         fdyn->mlvar->transterm=1;
      else if (strncmp(buffer,"neglect_lst",11)==0)
         fdyn->mlvar->transterm=2;
      else if (strncmp(buffer,"neglect_both",12)==0)
         fdyn->mlvar->transterm=3;
      else
         dserror("unknown treatment of transient term");
   }
   frchar("QUASTATBUB"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no",2)==0)
         fdyn->mlvar->quastabub=0;
      else if (strncmp(buffer,"yes",3)==0)
         fdyn->mlvar->quastabub=1;
      else
         dserror("unknown quasi-static bubbles");
   }
   frchar("CONVECVEL"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no_approx",9)==0)
         fdyn->mlvar->convel=0;
      else if (strncmp(buffer,"ls_approx",9)==0)
         fdyn->mlvar->convel=1;
      else
         dserror("unknown convective velocity");
   }
   frchar("SMSGVISC"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no",2)==0)
         fdyn->mlvar->smsgvi=0;
      else if (strncmp(buffer,"artificial",10)==0)
         fdyn->mlvar->smsgvi=1;
      else if (strncmp(buffer,"Smagorinsky",11)==0)
         fdyn->mlvar->smsgvi=2;
      else if (strncmp(buffer,"dynamic",7)==0)
         fdyn->mlvar->smsgvi=3;
      else
         dserror("Unknown submesh subgrid viscosity!");
   }
   frchar("SMSTABI"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no",2)==0)
         fdyn->mlvar->smstabi=0;
      else if (strncmp(buffer,"yes",3)==0)
         fdyn->mlvar->smstabi=1;
      else
         dserror("unknown submesh stabilization");
   }
   frchar("SMSTABDO"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"USFEM_instat",12)==0)
         fdyn->mlvar->smstado=-1;
      else if (strncmp(buffer,"USFEM_CDR",9)==0)
         fdyn->mlvar->smstado=-2;
      else if (strncmp(buffer,"USFEM_stat",10)==0)
         fdyn->mlvar->smstado=-3;
      else if (strncmp(buffer,"GLS-_instat",11)==0)
         fdyn->mlvar->smstado=1;
      else if (strncmp(buffer,"GLS-_CDR",8)==0)
         fdyn->mlvar->smstado=2;
      else if (strncmp(buffer,"GLS-_stat",9)==0)
         fdyn->mlvar->smstado=3;
      else
         dserror("unknown differential operator for submesh stabilization");
   }
   frchar("SMSTABNORM"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"L_1",3)==0)
         fdyn->mlvar->smstano=1;
      else if (strncmp(buffer,"L_2",3)==0)
         fdyn->mlvar->smstano=2;
      else
         dserror("unknown submesh stabilization norm");
   }
   frchar("SMSTABMK"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"-1",2)==0)
         fdyn->mlvar->smstamk=-1;
      else if (strncmp(buffer,"0",1)==0)
         fdyn->mlvar->smstamk=0;
      else if (strncmp(buffer,"1",1)==0)
         fdyn->mlvar->smstamk=1;
      else
         dserror("unknown submesh mk");
   }
   frchar("SMSTABNINTHS"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"at_center",9)==0)
         fdyn->mlvar->smstani=0;
      else if (strncmp(buffer,"every_intpt",11)==0)
         fdyn->mlvar->smstani=1;
      else
         dserror("unknown submesh number of int. points for elesize");
   }
   frchar("SMORDER"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"linear",6)==0)
         fdyn->mlvar->smorder=1;
      else if (strncmp(buffer,"quadratic",9)==0)
         fdyn->mlvar->smorder=2;
      else
         dserror("unknown submesh element order");
   }
   frchar("SMNONUNIF"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no",2)==0)
         fdyn->mlvar->smnunif=0;
      else if (strncmp(buffer,"shishkin_type",13)==0)
         fdyn->mlvar->smnunif=1;
      else if (strncmp(buffer,"quadratic",9)==0)
         fdyn->mlvar->smnunif=2;
      else if (strncmp(buffer,"cubic",5)==0)
         fdyn->mlvar->smnunif=3;
      else
         dserror("unknown sub-submesh element order");
   }
   frchar("SSMORD"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"linear",6)==0)
         fdyn->mlvar->ssmorder=1;
      else if (strncmp(buffer,"quadratic",9)==0)
         fdyn->mlvar->ssmorder=2;
      else
         dserror("unknown sub-submesh element order");
   }
   frchar("SSMNNUN"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (strncmp(buffer,"no",2)==0)
         fdyn->mlvar->ssmnunif=0;
      else if (strncmp(buffer,"shishkin_type",13)==0)
         fdyn->mlvar->ssmnunif=1;
      else if (strncmp(buffer,"quadratic",9)==0)
         fdyn->mlvar->ssmnunif=2;
      else if (strncmp(buffer,"cubic",5)==0)
         fdyn->mlvar->ssmnunif=3;
      else
         dserror("unknown sub-submesh element order");
   }
/*--------------read int */
   frint("SMELESIZE"      ,&(fdyn->mlvar->smesize)  ,&ierr);
   frint("SMSTABPAR"      ,&(fdyn->mlvar->smstapa)  ,&ierr);
   frint("SMELEMENTS"     ,&(fdyn->mlvar->smelenum) ,&ierr);
   frint("SMNUMGP"        ,&(fdyn->mlvar->smnumgp)  ,&ierr);
   frint("SSMELEM"        ,&(fdyn->mlvar->ssmelenum),&ierr);
   frint("SSMNGP"         ,&(fdyn->mlvar->ssmnumgp) ,&ierr);
/*--------------read double */
   frdouble("SMSMAGCON"   ,&(fdyn->mlvar->smsmagcon),&ierr);

   frread();
}
frrewind();
/*----------------------------------------------------------------------*/

end:
/*-------------------------------------------------- plausibility check */
#ifndef D_FSI
if (fdyn->freesurf>0)
  dserror("Freesurf requires FSI functions to compile in!\n");
#endif

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inpctr_dyn_fluid */

#endif

#ifdef D_FSI

/*!---------------------------------------------------------------------
\brief input of the FSI DYNAMIC block in the input-file

<pre>                                                         genk 09/02

In this routine the data in the FSI DYNAMIC block of the input file
are read and stored in fsidyn

</pre>
\param  *fsidyn 	  FSI_DATA       (o)
\return void

------------------------------------------------------------------------*/
void inpctr_dyn_fsi(FSI_DYNAMIC *fsidyn)
{
INT    ierr;
char   buffer[50];

#ifdef DEBUG
dstrc_enter("inpctr_dyn_fsi");
#endif

/*-------------------------------------------------- set default values */
fsidyn->ifsi=fsi_basic_sequ_stagg;             /* default: sequ. staggered */
fsidyn->ipre=1;             /* default: d(n) */
fsidyn->inrmfsi=1;
fsidyn->ichecke=0;          /* default: no energy check */
fsidyn->isdmax=1;
fsidyn->nstep=1;
fsidyn->itemax=1;
fsidyn->iale=1;             /* default pseudo structure */
fsidyn->uppss=1;
fsidyn->upres=1;
fsidyn->uprestart=1;
fsidyn->dt=ONE/TEN;
fsidyn->maxtime=1000.0;
fsidyn->entol=EPS6;
fsidyn->relax=ONE;
fsidyn->convtol=EPS6;
fsidyn->coupmethod=1;
fsidyn->coupforce = cf_stress; /* couple by stresses */

if (frfind("-FSI DYNAMIC")==0) goto end;
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("COUPALGO"  ,buffer    ,&ierr);

   if (ierr==1)
   {
      if (frwordcmp(buffer,"basic_sequ_stagg")==0)
         fsidyn->ifsi=fsi_basic_sequ_stagg;
      else if (frwordcmp(buffer,"sequ_stagg_pred")==0)
         fsidyn->ifsi=fsi_sequ_stagg_pred;
      else if (frwordcmp(buffer,"sequ_stagg_shift")==0)
         fsidyn->ifsi=fsi_sequ_stagg_shift;
      else if (frwordcmp(buffer,"iter_stagg_fixed_rel_param")==0)
         fsidyn->ifsi=fsi_iter_stagg_fixed_rel_param;
      else if (frwordcmp(buffer,"iter_stagg_AITKEN_rel_param")==0)
         fsidyn->ifsi=fsi_iter_stagg_AITKEN_rel_param;
      else if (frwordcmp(buffer,"iter_stagg_steep_desc")==0)
         fsidyn->ifsi=fsi_iter_stagg_steep_desc;
      else if (frwordcmp(buffer,"iter_stagg_CHEB_rel_param")==0)
         fsidyn->ifsi=fsi_iter_stagg_CHEB_rel_param;
      else if (frwordcmp(buffer,"iter_stagg_AITKEN_rel_force")==0)
         fsidyn->ifsi=fsi_iter_stagg_AITKEN_rel_force;
      else if (frwordcmp(buffer,"iter_stagg_steep_desc_force")==0)
         fsidyn->ifsi=fsi_iter_stagg_steep_desc_force;
      else if (frwordcmp(buffer,"iter_stagg_Newton_FD")==0)
        fsidyn->ifsi=fsi_iter_stagg_Newton_FD;
      else if (frwordcmp(buffer,"iter_stagg_Newton_I")==0)
        fsidyn->ifsi=fsi_iter_stagg_Newton_I;
      else
         dserror("Coupling Algorithm COUPALGO unknown");
   }

   frchar("PREDICTOR"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"d(n)")==0)
         fsidyn->ipre=1;
      else if (frwordcmp(buffer,"d(n)+dt*(1.5*v(n)-0.5*v(n-1))")==0)
         fsidyn->ipre=2;
      else if (frwordcmp(buffer,"d(n)+dt*v(n)")==0)
            fsidyn->ipre=3;
      else if (frwordcmp(buffer,"d(n)+dt*v(n)+0.5*dt^2*a(n)")==0)
            fsidyn->ipre=4;
      else dserror("PREDICTOR unknown!\n");
   }
   frchar("CONVCRIT"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"||g(i)||:sqrt(neq)")==0)
         fsidyn->inrmfsi=1;
      else if (frwordcmp(buffer,"||g(i)||:||g(0)||")==0)
         fsidyn->inrmfsi=2;
      else
         dserror("Convergence criterion CONVCRIT unknown!");
   }
   frreadyes("ENERGYCHECK", &(fsidyn->ichecke));

   frchar("IALE"    ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"Pseudo_Structure")==0)
         fsidyn->iale=1;
      else
         dserror("Parameter unknown: IALE");
   }
   /* parameter to indicate the coupling method that is applied         */
   /* mortar method       --> fsidyn-->coupmethod == 0                  */
   /* conforming coupling --> fsidyn-->coupmethod == 1 (default)        */
   frchar("COUPMETHOD"    ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"MTR")==0 ||
          frwordcmp(buffer,"Mtr")==0 ||
          frwordcmp(buffer,"mtr")==0)
          fsidyn->coupmethod=0;
      else if (frwordcmp(buffer,"conforming")==0)
               fsidyn->coupmethod=1;
      else
         dserror("Parameter unknown: COUPMETHOD");
   }
   frchar("COUPFORCE"    ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"none")==0)
          fsidyn->coupforce=cf_none;
      else if (frwordcmp(buffer,"stress")==0)
          fsidyn->coupforce=cf_stress;
      else if (frwordcmp(buffer,"nodeforce")==0)
          fsidyn->coupforce=cf_nodeforce;
      else
         dserror("Parameter unknown: COUPFORCE");
   }
/*--------------read INT */
   frint("ISDMAX"     ,&(fsidyn->isdmax)     ,&ierr);
   frint("NUMSTEP"    ,&(fsidyn->nstep)      ,&ierr);
   frint("ITEMAX"     ,&(fsidyn->itemax)     ,&ierr);
   frint("UPPSS"      ,&(fsidyn->uppss)      ,&ierr);
   frint("UPRES"      ,&(fsidyn->upres)      ,&ierr);
   frint("RESTARTEVRY",&(fsidyn->uprestart)  ,&ierr);
/*--------------read DOUBLE */
   frdouble("TIMESTEP"    ,&(fsidyn->dt)     ,&ierr);
   frdouble("MAXTIME"     ,&(fsidyn->maxtime),&ierr);
   frdouble("TOLENCHECK"  ,&(fsidyn->entol)  ,&ierr);
   frdouble("RELAX"       ,&(fsidyn->relax)  ,&ierr);
   frdouble("CONVTOL"     ,&(fsidyn->convtol),&ierr);
   frread();
}
frrewind();

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inpctr_dyn_fsi */

#endif

#ifdef D_SSI

/*!---------------------------------------------------------------------
\brief input of the SSI DYNAMIC block in the input-file

<pre>                                                         genk 10/03

In this routine the data in the SSI DYNAMIC block of the input file
are read and stored in ssidyn

</pre>
\param  *ssidyn 	  SSI_DATA       (o)
\return void

------------------------------------------------------------------------*/
void inpctr_dyn_ssi(SSI_DYNAMIC *ssidyn)
{
INT    ierr;
INT    mod_res_write;
char   buffer[50];
#ifdef DEBUG
dstrc_enter("inpctr_dyn_ssi");
#endif

/*-------------------------------------------------- set default values */
ssidyn->ifsi=1;
ssidyn->ipre=1;
ssidyn->inrmfsi=1;
ssidyn->ichecke=0;
ssidyn->inest=0;
ssidyn->ichopt=0;
ssidyn->iait=0;
ssidyn->itechapp=0;
ssidyn->ichmax=1;
ssidyn->isdmax=1;
ssidyn->nstep=1;
ssidyn->itemax=1;
ssidyn->iale=1;
ssidyn->uppss=1;
ssidyn->upres=1;
ssidyn->res_write_evry=1;
ssidyn->dt=ONE/TEN;
ssidyn->maxtime=1000.0;
ssidyn->entol=EPS6;
ssidyn->relax=ONE;
ssidyn->convtol=EPS6;
ssidyn->conformmesh=1;
ssidyn->coupmethod=0;

if (frfind("-SSI DYNAMIC")==0) goto end;
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("COUPALGO"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"basic_sequ_stagg")==0)
         ssidyn->ifsi=1;
      else if (frwordcmp(buffer,"sequ_stagg_pred")==0)
         ssidyn->ifsi=2;
      else if (frwordcmp(buffer,"sequ_stagg_shift")==0)
         ssidyn->ifsi=3;
      else if (frwordcmp(buffer,"iter_stagg_fixed_rel_param")==0)
         ssidyn->ifsi=4;
      else if (frwordcmp(buffer,"iter_stagg_AITKEN_rel_param")==0)
         ssidyn->ifsi=5;
      else if (frwordcmp(buffer,"iter_stagg_steep_desc")==0)
         ssidyn->ifsi=6;
      else if (frwordcmp(buffer,"iter_stagg_CHEB_rel_param")==0)
         ssidyn->ifsi=7;
      else
         dserror("Coupling Algorithm COUPALGO unknown");
   }
   frchar("CONVCRIT"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"||g(i)||:sqrt(neq)")==0)
         ssidyn->inrmfsi=1;
      else if (frwordcmp(buffer,"||g(i)||:||g(0)||")==0)
         ssidyn->inrmfsi=2;
      else
         dserror("Convergence criterion CONVCRIT unknown!");
   }
   frreadyes("ENERGYCHECK",&(ssidyn->ichecke));

   frchar("NESTEDITER"  ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frcheckno(buffer))
         ssidyn->inest=0;
      else if (frwordcmp(buffer,"fluid+structure")==0)
         ssidyn->inest=1;
      else if (frwordcmp(buffer,"structure")==0)
         ssidyn->inest=2;
      else
         dserror("Parameter unknown: NESTEDITER");
   }
   frreadyes("ICHOPT",&(ssidyn->ichopt));
   frchar("IALE"    ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"Pseudo_Structure")==0)
         ssidyn->iale=1;
      else
         dserror("Parameter unknown: IALE");
   }

   /* This flag is used the wrong way. Just to keep the code
    * interesting. */
   if (frreadyes("MESHCONFORM",&(ssidyn->conformmesh)))
   {
     ssidyn->conformmesh = !ssidyn->conformmesh;
   }

   /* parameter to indicate the coupling method that is applied */
   /* mortar method --> ssidyn->coupmethod == 0 (default) */
   /* interpolation scheme --> ssidyn-->coupmethod == 1*/
   frchar("COUPMETHOD"    ,buffer    ,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"MTR")==0 ||
          frwordcmp(buffer,"Mtr")==0 ||
          frwordcmp(buffer,"mtr")==0)
          ssidyn->coupmethod=0;
      else if (frwordcmp(buffer,"INT")==0 ||
               frwordcmp(buffer,"Int")==0 ||
               frwordcmp(buffer,"int")==0)
               ssidyn->coupmethod=1;
      else
         dserror("Parameter unknown: COUPMETHOD");
   }
/*--------------read INT */
   frint("ITECHAPP"   ,&(ssidyn->itechapp)         ,&ierr);
   frint("ICHMAX"     ,&(ssidyn->ichmax)           ,&ierr);
   frint("ISDMAX"     ,&(ssidyn->isdmax)           ,&ierr);
   frint("NUMSTEP"    ,&(ssidyn->nstep)            ,&ierr);
   frint("ITEMAX"     ,&(ssidyn->itemax)           ,&ierr);
   frint("UPPSS"      ,&(ssidyn->uppss)            ,&ierr);
   frint("UPRES"      ,&(ssidyn->upres)            ,&ierr);
   frint("RESTARTEVRY",&(ssidyn->res_write_evry)   ,&ierr);
/*--------------read DOUBLE */
   frdouble("TIMESTEP"   ,&(ssidyn->dt)     ,&ierr);
   frdouble("MAXTIME"    ,&(ssidyn->maxtime),&ierr);
   frdouble("TOLENCHECK" ,&(ssidyn->entol)  ,&ierr);
   frdouble("RELAX"      ,&(ssidyn->relax)  ,&ierr);
   frdouble("CONVTOL"    ,&(ssidyn->convtol),&ierr);
   frread();
}
frrewind();
/*----------------------------------------------------------------------*/

end:
/*------------------------------------------------ check restart option */
if (ssidyn->res_write_evry >= ssidyn->upres)
{
   mod_res_write = ssidyn->res_write_evry % ssidyn->upres;
   if (mod_res_write!=0)
   dserror("RESTARTEVRY and UPRES have to be multiple of each other!\n");
}
else
{
   mod_res_write =  ssidyn->upres % ssidyn->res_write_evry;
   if (mod_res_write!=0)
   dserror("RESTARTEVRY and UPRES have to be multiple of each other!\n");
}
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inpctr_dyn_ssi */

#endif


/*!---------------------------------------------------------------------
\brief input of the ALE DYNAMIC block in the input-file

<pre>                                                           ck 12/02

In this routine the data in the ALE DYNAMIC block of the input file
are read and stored in adyn

</pre>
\param  *adyn 	  ALE_DATA       (o)
\return void

------------------------------------------------------------------------*/
#ifdef D_ALE
void inpctr_dyn_ale(ALE_DYNAMIC *adyn)
{
INT    ierr;
char buffer[50];

#ifdef DEBUG
dstrc_enter("inpctr_dyn_ale");
#endif
/*---------------------------------------------------- some defaults ---*/
adyn->num_initstep = 0;
adyn->step = 0;
adyn->nstep = 10;
adyn->dt = 0.1;
adyn->maxtime = 1000.0;
adyn->updevry_disp = 20;
adyn->measure_quality = no_quality;
adyn->typ = classic_lin;

if (frfind("-ALE DYNAMIC")==0) goto end;
frread();
while(strncmp(allfiles.actplace,"------",6)!=0)
{
/*--------------read chars */
   frchar("ALE_TYPE",buffer,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"classic_lin")==0) adyn->typ = classic_lin;
      else if (frwordcmp(buffer,"min_Je_stiff")==0) adyn->typ = min_Je_stiff;
      else if (frwordcmp(buffer,"two_step")==0) adyn->typ = two_step;
      else if (frwordcmp(buffer,"springs")==0) adyn->typ = springs;
      else if (frwordcmp(buffer,"laplace")==0) adyn->typ = laplace;
      else if (frwordcmp(buffer,"LAS")==0) adyn->typ = LAS;
      else dserror("unknown ALE_TYPE");
   }
   frchar("QUALITY",buffer,&ierr);
   if (ierr==1)
   {
      if (frwordcmp(buffer,"aspect_ratio")==0) adyn->measure_quality = aspect_ratio;
      else if (frwordcmp(buffer,"ASPECT_RATIO")==0) adyn->measure_quality = aspect_ratio;
      else if (frwordcmp(buffer,"corner_angle")==0) adyn->measure_quality = corner_angle;
      else if (frwordcmp(buffer,"CORNER_ANGLE")==0) adyn->measure_quality = corner_angle;
      else if (frwordcmp(buffer,"min_Je")==0) adyn->measure_quality = min_detF;
      else if (frwordcmp(buffer,"MIN_Je")==0) adyn->measure_quality = min_detF;
      else if (frwordcmp(buffer,"none")==0) adyn->measure_quality = no_quality;
      else if (frwordcmp(buffer,"NONE")==0) adyn->measure_quality = no_quality;
      else dserror("unknown QUALITY in ALE DYNAMIC");
   }
/*--------------read INT */
   frint("NUMSTEP" ,   &(adyn->nstep) ,       &ierr);
   frint("RESEVRYDISP",&(adyn->updevry_disp), &ierr);
   frint("NUM_INITSTEP",&(adyn->num_initstep), &ierr);
/*--------------read DOUBLE */
   frdouble("TIMESTEP", &(adyn->dt)     , &ierr);
   frdouble("MAXTIME" , &(adyn->maxtime), &ierr);

   frread();
}
frrewind();

end:
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inpctr_dyn_ale */
#endif




/*----------------------------------------------------------------------*/
/*!
\brief Input of the TSI DYNAMIC block in the input-file

In this routine the data in the TSI DYNAMIC block of the input file
are read and stored in tsidyn

\param  *tsidyn 	  TSI_DYNAMIC       (o)
\return void

\author bborn
\date 03/06
*/
/*------------------------------------------------------------------------*/
#ifdef D_TSI
void inpctr_dyn_tsi(TSI_DYNAMIC *tsidyn)
{
  INT    ierr;
  char   buffer[50];

#ifdef DEBUG
  dstrc_enter("inpctr_dyn_tsi");
#endif

  /*--------------------------------------------------------------------*/
  /* set important defaults */
  tsidyn->kind = tsi_full;  /* default: full thermostructural coupling */

  /*--------------------------------------------------------------------*/
  /* read settings in input file */
  if (frfind("-TSI DYNAMIC") == 0) goto end;
  frread();
  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    /* read chars */
    frchar("KIND", buffer, &ierr);
    if (ierr == 1)
    {
      if (frwordcmp(buffer,"full") == 0)
      {
	 tsidyn->kind = tsi_full;
      }
      else if (frwordcmp(buffer, "thermal_static_structure_dynamic") == 0)
      {
	 tsidyn->kind = tsi_therm_stat_struct_dyn;
      }
      else if (frwordcmp(buffer, "thermal_predefined_structure_dynamic") == 0)
      {
         tsidyn->kind = tsi_therm_pred_struct_dyn;
      }
      else
      {
         dserror("KIND of TSI unknown");
      }
    }

    frint("NUMSTEP", &(tsidyn->nstep), &ierr);
    frint("TIMESTEP", &(tsidyn->dt), &ierr);
    frdouble("MAXTIME", &(tsidyn->maxtime), &ierr);

    frint("ITEMAX", &(tsidyn->maxiter), &ierr);

    frread();
  }
  frrewind();

  end:
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of inpctr_dyn_tsi */
#endif  /* end of #ifdef D_TSI */
